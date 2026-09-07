#!/usr/bin/env python3

import argparse
import re
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.patches as patches


DEFAULT_MOTIFS = "AAAGGG,AAAAGG,AAAAG,AAAGG,AAGGG,ACAGG,AAAGC"

MOTIF_COLORS = {
    "AAAGGG": "#0072B2",
    "AAAAGG": "#E69F00",
    "AAAAG": "#009E73",
    "AAAGG": "#D55E00",
    "AAGGG": "#CC79A7",
    "ACAGG": "#8C564B",
    "AAAGC": "#E377C2",
    "OTHER": "#BDBDBD",
}


def read_fastq(path):
    with open(path) as f:
        while True:
            header = f.readline().rstrip("\n")
            if not header:
                break
            seq = f.readline().rstrip("\n").upper()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                break
            yield header, seq


def parse_header(header):
    h = header[1:] if header.startswith("@") else header
    first = h.split()[0]

    if first.startswith("STD_"):
        read_type = "STD"
        label = first.split("|")[0]
    elif first.startswith("SC_L_"):
        read_type = "SC_L"
        label = first.split("|")[0]
    elif first.startswith("SC_R_"):
        read_type = "SC_R"
        label = first.split("|")[0]
    elif first.startswith("INS_"):
        read_type = "INS"
        label = first.split("|")[0]
    else:
        read_type = "READ"
        label = first.split("|")[0]

    flag = 0
    m = re.search(r"flag=(\d+)", h)
    if m:
        flag = int(m.group(1))

    ref_pos = None
    m = re.search(r"pos=([^|]+)", h)
    if m:
        ref_pos = m.group(1)

    return {
        "label": label,
        "type": read_type,
        "flag": flag,
        "ref_pos": ref_pos,
        "header": header,
    }


def ref_pos_to_int(ref_pos):
    if ref_pos is None:
        return None
    m = re.search(r":(\d+)", ref_pos)
    if not m:
        return None
    return int(m.group(1))


def infer_plot_side(info, repeat_start, repeat_end):
    """
    Infer plotting side for two-sided display.

    If reference position is available:
      positions before the repeat core are placed on the left side;
      positions after the repeat core are placed on the right side.

    If the position falls within the repeat core or is unavailable,
    fall back to the read label.
    """
    pos = ref_pos_to_int(info.get("ref_pos"))

    if pos is not None:
        if pos < repeat_start:
            return "L"
        if pos > repeat_end:
            return "R"

    if info["type"] == "SC_L":
        return "R"
    if info["type"] == "SC_R":
        return "L"

    return "L"


def scan_motifs(seq, motifs):
    """
    Greedy motif scan from left to right.
    Non-matching bases are grouped as OTHER.
    """
    i = 0
    blocks = []

    while i < len(seq):
        matched = None
        for motif in motifs:
            if seq.startswith(motif, i):
                matched = motif
                break

        if matched:
            start = i
            i += len(matched)
            while seq.startswith(matched, i):
                i += len(matched)
            blocks.append((start, i, matched))
        else:
            start = i
            i += 1
            while i < len(seq):
                if any(seq.startswith(m, i) for m in motifs):
                    break
                i += 1
            blocks.append((start, i, "OTHER"))

    return blocks


def motif_summary(blocks):
    counts = Counter()
    for start, end, motif in blocks:
        counts[motif] += end - start

    if counts:
        dominant_motif, dominant_bp = counts.most_common(1)[0]
    else:
        dominant_motif, dominant_bp = "NA", 0

    return counts, dominant_motif, dominant_bp


def short_label(label, max_len=18):
    if len(label) <= max_len:
        return label
    return label[:max_len]


def main():
    parser = argparse.ArgumentParser(
        description="Generate RFC1 repeat motif composition plots from FASTQ reads."
    )

    parser.add_argument("--fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--out", required=True, help="Output prefix")
    parser.add_argument("--motifs", default=DEFAULT_MOTIFS,
                        help="Comma-separated repeat motifs")
    parser.add_argument("--repeat-start", type=int, default=39348425,
                        help="RFC1 repeat core start, 1-based inclusive")
    parser.add_argument("--repeat-end", type=int, default=39348479,
                        help="RFC1 repeat core end, 1-based inclusive")
    parser.add_argument("--title", default="RFC1 repeat motif composition plot",
                        help="Plot title")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Output image resolution")
    parser.add_argument("--width", type=float, default=12,
                        help="Figure width")
    parser.add_argument("--row-height", type=float, default=0.35,
                        help="Height per read row")
    parser.add_argument("--min-height", type=float, default=4.0,
                        help="Minimum figure height")

    args = parser.parse_args()

    motifs = [m.strip().upper() for m in args.motifs.split(",") if m.strip()]
    motifs = sorted(motifs, key=len, reverse=True)

    records = []

    for header, seq in read_fastq(args.fastq):
        info = parse_header(header)
        blocks = scan_motifs(seq, motifs)
        counts, dominant_motif, dominant_bp = motif_summary(blocks)
        side = infer_plot_side(info, args.repeat_start, args.repeat_end)

        records.append({
            "label": info["label"],
            "type": info["type"],
            "side": side,
            "flag": info["flag"],
            "length": len(seq),
            "blocks": blocks,
            "dominant_motif": dominant_motif,
            "dominant_bp": dominant_bp,
            "counts": counts,
            "header": header,
        })

    if not records:
        raise SystemExit("No reads found in input FASTQ.")

    # Sort STD reads first, then by side and length.
    type_order = {"STD": 0, "INS": 1, "SC_L": 2, "SC_R": 2, "READ": 3}
    records.sort(
        key=lambda r: (
            type_order.get(r["type"], 9),
            r["side"],
            -r["length"],
            r["label"]
        )
    )

    max_len = max(r["length"] for r in records)

    fig_height = max(args.min_height, len(records) * args.row_height + 2.2)
    fig, ax = plt.subplots(figsize=(args.width, fig_height))

    placement_rows = []

    for y, rec in enumerate(records):
        if rec["side"] == "R":
            x_offset = max_len - rec["length"]
        else:
            x_offset = 0

        arrow = "←" if rec["flag"] & 16 else "→"
        y_plot = len(records) - 1 - y

        for start, end, motif in rec["blocks"]:
            x = x_offset + start
            width = end - start
            color = MOTIF_COLORS.get(motif, MOTIF_COLORS["OTHER"])

            rect = patches.Rectangle(
                (x, y_plot - 0.4),
                width,
                0.8,
                linewidth=0,
                facecolor=color,
            )
            ax.add_patch(rect)

        label = f"{short_label(rec['label'])} {arrow}"
        ax.text(
            -max_len * 0.015,
            y_plot,
            label,
            va="center",
            ha="right",
            fontsize=7,
        )

        placement_rows.append({
            "label": f"{rec['label']} {arrow}",
            "type": rec["type"],
            "side": rec["side"],
            "flag": rec["flag"],
            "length_bp": rec["length"],
            "x_start": x_offset,
            "x_end": x_offset + rec["length"],
            "dominant_motif": rec["dominant_motif"],
            "dominant_bp": rec["dominant_bp"],
            "AAGGG_run_bp": rec["counts"].get("AAGGG", 0),
            "AAAGG_run_bp": rec["counts"].get("AAAGG", 0),
            "original_header": rec["header"],
        })

    ax.set_xlim(-max_len * 0.08, max_len)
    ax.set_ylim(-0.8, len(records) - 0.2)

    ax.set_xlabel("Position from the left reference-side anchor (bp)")
    ax.set_ylabel("Read ID")
    ax.set_title(args.title)

    ax.set_yticks([])

    ax.text(0, len(records) - 0.05, "Reference-left side",
            ha="left", va="bottom", fontsize=8)
    ax.text(max_len, len(records) - 0.05, "Reference-right side",
            ha="right", va="bottom", fontsize=8)

    legend_order = [
        "AAAGGG", "AAAAGG", "AAAAG", "AAAGG",
        "AAGGG", "ACAGG", "AAAGC", "OTHER"
    ]
    handles = []
    labels = []

    for motif in legend_order:
        if motif in MOTIF_COLORS:
            handles.append(
                patches.Patch(color=MOTIF_COLORS[motif], label=motif)
            )
            labels.append(motif)

    ax.legend(
        handles=handles,
        labels=labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=4,
        fontsize=7,
        frameon=False,
    )

    plt.tight_layout()
    png_path = f"{args.out}.png"
    tsv_path = f"{args.out}.placement.tsv"

    fig.savefig(png_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    with open(tsv_path, "w") as out:
        out.write(
            "label\ttype\tside\tflag\tlength_bp\tx_start\tx_end\t"
            "dominant_motif\tdominant_bp\tAAGGG_run_bp\tAAAGG_run_bp\toriginal_header\n"
        )
        for row in placement_rows:
            out.write(
                f"{row['label']}\t{row['type']}\t{row['side']}\t"
                f"{row['flag']}\t{row['length_bp']}\t{row['x_start']}\t"
                f"{row['x_end']}\t{row['dominant_motif']}\t"
                f"{row['dominant_bp']}\t{row['AAGGG_run_bp']}\t"
                f"{row['AAAGG_run_bp']}\t{row['original_header']}\n"
            )

    print(f"PNG: {png_path}")
    print(f"TSV: {tsv_path}")
    print(f"Reads: {len(records)}")


if __name__ == "__main__":
    main()
