#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import yaml


QUERY_LINE_RE = re.compile(r"^Query\s+(?P<start>\d+)\s+(?P<seq>[A-Za-z\-]+)\s+(?P<end>\d+)\s*$")
SBJCT_LINE_RE = re.compile(r"^Sbjct\s+(?P<start>\d+)\s+(?P<seq>[A-Za-z\-]+)\s+(?P<end>\d+)\s*$")
SEQ_ID_RE = re.compile(r"^Sequence ID:\s*(?P<seqid>\S+)\s+Length:\s*(?P<length>\S+)\s*$")
RANGE_RE = re.compile(r"^Range\s+(?P<range_no>\d+):\s*(?P<start>\S+)\s+to\s+(?P<end>\S+)\s*$")


SUMMARY_FIELD_ORDER = [
    "description",
    "scientific_name",
    "common_name",
    "taxid",
    "max_score",
    "total_score",
    "query_cover",
    "evalue",
    "percent_ident",
    "length",
    "accession",
]


NUMERIC_FIELDS = {"taxid", "length"}
FLOAT_FIELDS = {"max_score", "total_score", "percent_ident"}


class BlastParseError(RuntimeError):
    pass


def load_yaml(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    if not isinstance(config, dict):
        raise BlastParseError("YAML のルートは mapping(dict) である必要があります。")
    return config


def load_lines(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        return [line.rstrip("\n") for line in fh]


def find_line_index(lines: List[str], keyword: str, start: int = 0, exact: bool = False) -> int:
    for i in range(start, len(lines)):
        text = lines[i].strip()
        if exact:
            if text == keyword:
                return i
        else:
            if keyword in lines[i]:
                return i
    raise BlastParseError(f"キーワードが見つかりません: {keyword!r}")


def previous_nonempty_index(lines: List[str], idx: int) -> int:
    for i in range(idx - 1, -1, -1):
        if lines[i].strip():
            return i
    raise BlastParseError("前方に非空行が見つかりません。")


def build_summary_slices(abbrev_header: str, detail_header: str) -> List[Tuple[str, int, Optional[int]]]:
    starts = {
        "description": detail_header.index("Description"),
        "scientific_name": abbrev_header.index("Scientific"),
        "common_name": abbrev_header.index("Common"),
        "taxid": detail_header.index("Taxid"),
        "max_score": abbrev_header.index("Max"),
        "total_score": abbrev_header.index("Total"),
        "query_cover": abbrev_header.index("Query"),
        "evalue": detail_header.index("Value"),
        "percent_ident": detail_header.index("Ident"),
        "length": detail_header.index("Len"),
        "accession": detail_header.index("Accession"),
    }
    ordered = sorted(starts.items(), key=lambda x: x[1])
    slices: List[Tuple[str, int, Optional[int]]] = []
    for i, (name, start) in enumerate(ordered):
        end = ordered[i + 1][1] if i + 1 < len(ordered) else None
        slices.append((name, start, end))
    return slices


def cast_summary_value(field: str, value: str) -> Any:
    value = value.strip()
    if value == "":
        return ""
    if field in NUMERIC_FIELDS:
        return int(value) if value.isdigit() else value
    if field in FLOAT_FIELDS:
        try:
            return float(value)
        except ValueError:
            return value
    return value


def parse_summary_line_by_slices(line: str, slices: List[Tuple[str, int, Optional[int]]]) -> Dict[str, Any]:
    record: Dict[str, Any] = {}
    for field, start, end in slices:
        raw = line[start:end].strip() if end is not None else line[start:].strip()
        record[field] = cast_summary_value(field, raw)
    return record


def parse_summary_section(lines: List[str], cfg: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    start_keyword = cfg.get("start_keyword", "Description")
    end_keyword = cfg.get("end_keyword", "Alignments:")
    key_field = cfg.get("key_field", "accession")

    header_idx = find_line_index(lines, start_keyword)
    abbrev_idx = previous_nonempty_index(lines, header_idx)
    end_idx = find_line_index(lines, end_keyword, start=header_idx + 1, exact=True)

    slices = build_summary_slices(lines[abbrev_idx], lines[header_idx])
    summary_dict: Dict[str, Dict[str, Any]] = {}

    for i in range(header_idx + 1, end_idx):
        line = lines[i]
        if not line.strip():
            continue
        record = parse_summary_line_by_slices(line, slices)
        key = str(record.get(key_field, "")).strip()
        if not key:
            continue
        summary_dict[key] = record
    return summary_dict


def parse_sequence_id(line: str, prefix: str) -> Optional[Tuple[str, Any]]:
    if not line.strip().startswith(prefix):
        return None
    m = SEQ_ID_RE.match(line.strip())
    if not m:
        payload = line.split(prefix, 1)[1].strip()
        parts = payload.split()
        if not parts:
            return None
        return parts[0], ""
    seqid = m.group("seqid")
    length = m.group("length")
    return seqid, int(length) if length.isdigit() else length


def parse_range_line(line: str) -> Dict[str, Any]:
    m = RANGE_RE.match(line.strip())
    if not m:
        return {"range_raw": line.strip()}
    result: Dict[str, Any] = {
        "range_raw": line.strip(),
        "range_no": int(m.group("range_no")),
        "range_start": m.group("start"),
        "range_end": m.group("end"),
    }
    if result["range_start"].isdigit():
        result["range_start"] = int(result["range_start"])
    if result["range_end"].isdigit():
        result["range_end"] = int(result["range_end"])
    return result


def parse_alignment_triplet(
    query_line: str,
    match_line: str,
    sbjct_line: str,
    current_range: Optional[Dict[str, Any]],
    current_score: Optional[str],
    current_descriptions: List[str],
) -> Dict[str, Any]:
    q = QUERY_LINE_RE.match(query_line.strip())
    s = SBJCT_LINE_RE.match(sbjct_line.strip())
    if not q:
        raise BlastParseError(f"Query 行を解釈できません: {query_line!r}")
    if not s:
        raise BlastParseError(f"Sbjct 行を解釈できません: {sbjct_line!r}")

    record: Dict[str, Any] = {
        "query_line_raw": query_line.rstrip(),
        "match_line_raw": match_line.rstrip(),
        "sbjct_line_raw": sbjct_line.rstrip(),
        "query_start": int(q.group("start")),
        "query_seq": q.group("seq"),
        "query_end": int(q.group("end")),
        "sbjct_start": int(s.group("start")),
        "sbjct_seq": s.group("seq"),
        "sbjct_end": int(s.group("end")),
        "hit_descriptions": list(current_descriptions),
        "score_raw": current_score or "",
    }
    if current_range:
        record.update(current_range)
    return record


def parse_alignment_section(lines: List[str], cfg: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    start_keyword = cfg.get("start_keyword", "Alignments:")
    sequence_id_prefix = cfg.get("sequence_id_prefix", "Sequence ID:")
    query_prefix = cfg.get("query_prefix", "Query")
    sbjct_prefix = cfg.get("sbjct_prefix", "Sbjct")

    start_idx = find_line_index(lines, start_keyword, exact=True)
    alignment_dict: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    current_seq_ids: List[str] = []
    current_descriptions: List[str] = []
    current_range: Optional[Dict[str, Any]] = None
    current_score: Optional[str] = None
    seen_alignment_for_current_group = False

    i = start_idx + 1
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if not stripped:
            i += 1
            continue

        if line.startswith(">"):
            if seen_alignment_for_current_group:
                current_seq_ids = []
                current_descriptions = []
                current_range = None
                current_score = None
                seen_alignment_for_current_group = False
            current_descriptions.append(line[1:].strip())
            i += 1
            continue

        seq_id_info = parse_sequence_id(line, sequence_id_prefix)
        if seq_id_info is not None:
            seq_id, seq_length = seq_id_info
            if seq_id not in current_seq_ids:
                current_seq_ids.append(seq_id)
            i += 1
            continue

        if stripped.startswith("Range "):
            current_range = parse_range_line(stripped)
            seen_alignment_for_current_group = True
            i += 1
            continue

        if stripped.startswith("Score:"):
            current_score = stripped
            i += 1
            continue

        if stripped.startswith(query_prefix):
            if i + 2 >= len(lines):
                raise BlastParseError("Query 行の後ろに十分な行がありません。")
            query_line = lines[i]
            match_line = lines[i + 1]
            sbjct_line = lines[i + 2]
            if not sbjct_line.strip().startswith(sbjct_prefix):
                raise BlastParseError(
                    "Query 行の 2 行後が Sbjct 行ではありません: "
                    f"{sbjct_line!r}"
                )
            record = parse_alignment_triplet(
                query_line=query_line,
                match_line=match_line,
                sbjct_line=sbjct_line,
                current_range=current_range,
                current_score=current_score,
                current_descriptions=current_descriptions,
            )
            if not current_seq_ids:
                raise BlastParseError(
                    "Sequence ID が見つかる前に Query/Sbjct セットが現れました。"
                )
            for seq_id in current_seq_ids:
                item = dict(record)
                item["sequence_id"] = seq_id
                alignment_dict[seq_id].append(item)
            i += 3
            continue

        i += 1

    return dict(alignment_dict)


def merge_summary_and_alignment(
    summary_dict: Dict[str, Dict[str, Any]],
    alignment_dict: Dict[str, List[Dict[str, Any]]],
) -> Dict[str, Dict[str, Any]]:
    merged: Dict[str, Dict[str, Any]] = {}
    for accession, summary in summary_dict.items():
        item = dict(summary)
        item["DIFF"] = alignment_dict.get(accession, [])
        merged[accession] = item
    return merged


def project_for_text(
    merged_dict: Dict[str, Dict[str, Any]],
    pickup_fields: List[str],
    diff_fields: List[str],
) -> List[Dict[str, Any]]:
    result: List[Dict[str, Any]] = []
    for accession, record in merged_dict.items():
        projected: Dict[str, Any] = {}
        for field in pickup_fields:
            if field == "DIFF":
                projected["DIFF"] = [
                    {diff_field: diff_item.get(diff_field, "") for diff_field in diff_fields}
                    for diff_item in record.get("DIFF", [])
                ]
            else:
                projected[field] = record.get(field, "")
        result.append(projected)
    return result


def build_table_rows(
    merged_dict: Dict[str, Dict[str, Any]],
    pickup_fields: List[str],
    diff_fields: List[str],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    scalar_fields = [field for field in pickup_fields if field != "DIFF"]
    include_diff = "DIFF" in pickup_fields

    for accession, record in merged_dict.items():
        diffs = record.get("DIFF", [])
        if include_diff and diffs:
            for idx, diff in enumerate(diffs, start=1):
                row = {field: record.get(field, "") for field in scalar_fields}
                row["diff_index"] = idx
                for field in diff_fields:
                    row[field] = diff.get(field, "")
                rows.append(row)
        else:
            row = {field: record.get(field, "") for field in scalar_fields}
            if include_diff:
                row["diff_index"] = ""
                for field in diff_fields:
                    row[field] = ""
            rows.append(row)
    return rows


def tsv_from_rows(rows: List[Dict[str, Any]]) -> str:
    if not rows:
        return ""
    headers = list(rows[0].keys())
    lines = ["\t".join(headers)]
    for row in rows:
        values: List[str] = []
        for key in headers:
            value = row.get(key, "")
            if isinstance(value, list):
                value = json.dumps(value, ensure_ascii=False)
            values.append(str(value))
        lines.append("\t".join(values))
    return "\n".join(lines) + "\n"


def _extract_match_seq(match_line_raw: str, expected_len: int) -> str:
    stripped = match_line_raw.rstrip("\n")
    if stripped.strip() == "":
        return " " * expected_len
    core = stripped.strip()
    if len(core) < expected_len:
        return core.ljust(expected_len)
    return core[:expected_len]


def _format_alignment_triplet(diff: Dict[str, Any]) -> str:
    query_seq = str(diff.get("query_seq", ""))
    sbjct_seq = str(diff.get("sbjct_seq", ""))
    seq_len = max(len(query_seq), len(sbjct_seq))
    query_seq = query_seq.ljust(seq_len)
    sbjct_seq = sbjct_seq.ljust(seq_len)
    match_seq = _extract_match_seq(str(diff.get("match_line_raw", "")), seq_len)

    query_start = str(diff.get("query_start", ""))
    query_end = str(diff.get("query_end", ""))
    sbjct_start = str(diff.get("sbjct_start", ""))
    sbjct_end = str(diff.get("sbjct_end", ""))

    start_w = max(len(query_start), len(sbjct_start), 1)
    end_w = max(len(query_end), len(sbjct_end), 1)
    left_label_w = 5  # len("Sbjct")
    mid_indent = " " * (left_label_w + 2 + start_w + 2)

    query_line = f"Query  {query_start:>{start_w}}  {query_seq}  {query_end:>{end_w}}"
    match_line = f"{mid_indent}{match_seq}"
    sbjct_line = f"Sbjct  {sbjct_start:>{start_w}}  {sbjct_seq}  {sbjct_end:>{end_w}}"
    return "\n".join([query_line, match_line, sbjct_line])


def monitor_text_from_records(
    merged_dict: Dict[str, Dict[str, Any]],
    pickup_fields: List[str],
) -> str:
    blocks: List[str] = []
    scalar_fields = [field for field in pickup_fields if field != "DIFF"]

    for accession, record in merged_dict.items():
        lines = [f"[{accession}]"]
        for field in scalar_fields:
            lines.append(f"{field}: {record.get(field, '')}")

        diffs = record.get("DIFF", [])
        if diffs:
            for idx, diff in enumerate(diffs, start=1):
                lines.append(f"--- alignment {idx} ---")
                lines.append(_format_alignment_triplet(diff))
        else:
            lines.append("(alignment not found)")
        blocks.append("\n".join(lines))

    return "\n\n".join(blocks) + "\n"


def render_output(
    merged_dict: Dict[str, Dict[str, Any]],
    output_cfg: Dict[str, Any],
    display_cfg: Dict[str, Any],
) -> str:
    mode = str(output_cfg.get("mode", "table")).lower()
    pickup_fields = display_cfg.get("pickup_fields", SUMMARY_FIELD_ORDER)
    diff_fields = display_cfg.get("diff_fields", [])

    if mode == "table":
        rows = build_table_rows(merged_dict, pickup_fields, diff_fields)
        return tsv_from_rows(rows)
    if mode in {"text", "txt"}:
        projected = project_for_text(merged_dict, pickup_fields, diff_fields)
        return json.dumps(projected, ensure_ascii=False, indent=2) + "\n"
    if mode == "monitor":
        return monitor_text_from_records(merged_dict, pickup_fields)
    raise BlastParseError(f"未対応の output.mode です: {mode!r}")


def parse_blast_text(config: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    input_cfg = config.get("input", {})
    parse_cfg = config.get("parse", {})
    blast_txt = input_cfg.get("blast_txt")
    if not blast_txt:
        raise BlastParseError("input.blast_txt が設定されていません。")

    lines = load_lines(Path(blast_txt))
    summary_dict = parse_summary_section(lines, parse_cfg.get("summary_section", {}))
    alignment_dict = parse_alignment_section(lines, parse_cfg.get("alignment_section", {}))
    merged_dict = merge_summary_and_alignment(summary_dict, alignment_dict)
    return merged_dict


def main() -> None:
    parser = argparse.ArgumentParser(description="BLASTP テキスト出力パーサ")
    parser.add_argument("config", help="YAML 設定ファイル")
    parser.add_argument(
        "--dump-merged-json",
        help="統合済み dict を JSON で保存するパス",
        default=None,
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    config = load_yaml(config_path)
    merged_dict = parse_blast_text(config)

    if args.dump_merged_json:
        Path(args.dump_merged_json).write_text(
            json.dumps(merged_dict, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    output_cfg = config.get("output", {})
    display_cfg = config.get("display", {})
    rendered = render_output(merged_dict, output_cfg, display_cfg)

    output_path = output_cfg.get("path")
    if output_path:
        Path(output_path).write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
