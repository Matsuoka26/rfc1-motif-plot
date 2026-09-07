#!/usr/bin/env python3

import argparse
import pysam


def phred_to_fastq(qualities):
    if qualities is None:
        return None
    return "".join(chr(q + 33) for q in qualities)


def write_fastq(handle, header, seq, qual=None):
    if not seq:
        return
    if qual is None:
        qual = "I" * len(seq)
    handle.write(f"@{header}\n")
    handle.write(f"{seq}\n")
    handle.write("+\n")
    handle.write(f"{qual}\n")


def terminal_softclips(aln):
    left_clip = 0
    right_clip = 0

    if aln.cigartuples:
        if aln.cigartuples[0][0] == 4:
            left_clip = aln.cigartuples[0][1]
        if aln.cigartuples[-1][0] == 4:
            right_clip = aln.cigartuples[-1][1]

    return left_clip, right_clip


def extract_query_interval(aln, ref_start_1based, ref_end_1based):
    """
    Extract the query sequence corresponding to a 1-based inclusive reference interval.
    Insertions within the interval are included.
    """
    seq = aln.query_sequence
    quals = phred_to_fastq(aln.query_qualities)

    if seq is None or aln.cigartuples is None:
        return None, None

    qpos = 0
    rpos = aln.reference_start + 1  # 1-based

    q_start = None
    q_end = None

    for op, length in aln.cigartuples:
        # M, =, X
        if op in (0, 7, 8):
            block_ref_start = rpos
            block_ref_end = rpos + length - 1

            if block_ref_end >= ref_start_1based and block_ref_start <= ref_end_1based:
                ov_start = max(ref_start_1based, block_ref_start)
                ov_end = min(ref_end_1based, block_ref_end)

                local_q_start = qpos + (ov_start - block_ref_start)
                local_q_end = qpos + (ov_end - block_ref_start) + 1

                if q_start is None:
                    q_start = local_q_start
                q_end = local_q_end

            qpos += length
            rpos += length

        # insertion to reference
        elif op == 1:
            # Include insertion if it falls inside the target reference interval
            if ref_start_1based <= rpos <= ref_end_1based + 1:
                if q_start is None:
                    q_start = qpos
                q_end = qpos + length
            qpos += length

        # deletion or skipped region
        elif op in (2, 3):
            rpos += length

        # soft clip
        elif op == 4:
            qpos += length

        # hard clip / padding
        elif op in (5, 6):
            pass

    if q_start is None or q_end is None or q_end <= q_start:
        return None, None

    out_seq = seq[q_start:q_end]
    out_qual = quals[q_start:q_end] if quals is not None else None

    return out_seq, out_qual


def main():
    parser = argparse.ArgumentParser(
        description="Extract RFC1 repeat-containing spanning and soft-clipped reads from a BAM file."
    )

    parser.add_argument("--bam", required=True, help="Input sorted BAM file")
    parser.add_argument("--chrom", default="chr4", help="Chromosome name, default: chr4")

    parser.add_argument("--core-start", type=int, default=39348425,
                        help="RFC1 repeat core start, 1-based inclusive")
    parser.add_argument("--core-end", type=int, default=39348479,
                        help="RFC1 repeat core end, 1-based inclusive")

    parser.add_argument("--search-start", type=int, default=39348119,
                        help="Surrounding search region start, 1-based inclusive")
    parser.add_argument("--search-end", type=int, default=39348876,
                        help="Surrounding search region end, 1-based inclusive")

    parser.add_argument("--min-softclip", type=int, default=100,
                        help="Minimum terminal soft-clip length to output as SC read")
    parser.add_argument("--max-terminal-softclip-for-std", type=int, default=100,
                        help="Maximum terminal soft-clip length allowed for STD classification")

    parser.add_argument("--include-supplementary", action="store_true",
                        help="Include supplementary alignments")
    parser.add_argument("--include-secondary", action="store_true",
                        help="Include secondary alignments")

    parser.add_argument("--out-prefix", required=True,
                        help="Output prefix")

    args = parser.parse_args()

    out_all = f"{args.out_prefix}.repeat_reads.fastq"
    out_std = f"{args.out_prefix}.STD.fastq"
    out_sc = f"{args.out_prefix}.SC.fastq"
    out_summary = f"{args.out_prefix}.summary.tsv"

    n_total = 0
    n_std = 0
    n_sc = 0

    with pysam.AlignmentFile(args.bam, "rb") as bam, \
            open(out_all, "w") as fh_all, \
            open(out_std, "w") as fh_std, \
            open(out_sc, "w") as fh_sc, \
            open(out_summary, "w") as summary:

        summary.write(
            "read_id\ttype\tflag\tref_start\tref_end\tleft_softclip\tright_softclip\toutput_length\n"
        )

        for aln in bam.fetch(
            args.chrom,
            args.search_start - 1,
            args.search_end
        ):
            if aln.is_unmapped:
                continue
            if aln.is_secondary and not args.include_secondary:
                continue
            if aln.is_supplementary and not args.include_supplementary:
                continue
            if aln.query_sequence is None:
                continue

            n_total += 1

            ref_start = aln.reference_start + 1
            ref_end = aln.reference_end

            left_clip, right_clip = terminal_softclips(aln)

            aln_kind = "supplementary" if aln.is_supplementary else "primary"

            # Spanning read: alignment covers both sides of the repeat core
            is_spanning = (
                ref_start <= args.core_start and
                ref_end >= args.core_end and
                left_clip <= args.max_terminal_softclip_for_std and
                right_clip <= args.max_terminal_softclip_for_std
            )

            if is_spanning:
                seq, qual = extract_query_interval(aln, args.core_start, args.core_end)

                if seq:
                    header = (
                        f"STD_{aln.query_name}|{aln_kind}|flag={aln.flag}|"
                        f"pos={args.chrom}:{ref_start}|end={args.chrom}:{ref_end}"
                    )
                    write_fastq(fh_all, header, seq, qual)
                    write_fastq(fh_std, header, seq, qual)

                    summary.write(
                        f"{aln.query_name}\tSTD\t{aln.flag}\t{ref_start}\t{ref_end}\t"
                        f"{left_clip}\t{right_clip}\t{len(seq)}\n"
                    )
                    n_std += 1

            # Soft-clipped reads
            seq = aln.query_sequence
            qual = phred_to_fastq(aln.query_qualities)

            if left_clip >= args.min_softclip:
                sc_seq = seq[:left_clip]
                sc_qual = qual[:left_clip] if qual is not None else None

                header = (
                    f"SC_L_{aln.query_name}|{aln_kind}|flag={aln.flag}|clip={left_clip}|"
                    f"pos={args.chrom}:{ref_start}|end={args.chrom}:{ref_end}"
                )
                write_fastq(fh_all, header, sc_seq, sc_qual)
                write_fastq(fh_sc, header, sc_seq, sc_qual)

                summary.write(
                    f"{aln.query_name}\tSC_L\t{aln.flag}\t{ref_start}\t{ref_end}\t"
                    f"{left_clip}\t{right_clip}\t{len(sc_seq)}\n"
                )
                n_sc += 1

            if right_clip >= args.min_softclip:
                sc_seq = seq[-right_clip:]
                sc_qual = qual[-right_clip:] if qual is not None else None

                header = (
                    f"SC_R_{aln.query_name}|{aln_kind}|flag={aln.flag}|clip={right_clip}|"
                    f"pos={args.chrom}:{ref_start}|end={args.chrom}:{ref_end}"
                )
                write_fastq(fh_all, header, sc_seq, sc_qual)
                write_fastq(fh_sc, header, sc_seq, sc_qual)

                summary.write(
                    f"{aln.query_name}\tSC_R\t{aln.flag}\t{ref_start}\t{ref_end}\t"
                    f"{left_clip}\t{right_clip}\t{len(sc_seq)}\n"
                )
                n_sc += 1

    print("Done.")
    print(f"Input alignments checked: {n_total}")
    print(f"STD reads written: {n_std}")
    print(f"SC reads written: {n_sc}")
    print(f"Output FASTQ: {out_all}")
    print(f"STD FASTQ: {out_std}")
    print(f"SC FASTQ: {out_sc}")
    print(f"Summary: {out_summary}")


if __name__ == "__main__":
    main()
