import argparse
import gzip


def load_fastq(fq):
    """Load a fastq file and return a list of lines."""
    if fq.endswith(".gz"):
        with gzip.open(fq, "rt") as f:
            return [line.strip() for line in f]
    else:
        with open(fq, "r") as f:
            return [line.strip() for line in f]


def trim_adapter(loaded_fq, adapter, min_len, trim5, trim3, min_5match, error):
    """Trim the adapters from the reads, 5'3 and 3' clip if needed, and return a list of trimmed reads."""
    trimmed_reads = []
    adapter_missing_count = 0
    too_short_count = 0
    for read_no in range(1, len(loaded_fq), 4):
        read = loaded_fq[read_no]
        if adapter[:min_5match] in read:
            if read.index(adapter[:min_5match]) >= min_len + trim5 + trim3:
                trimmed_read = read[trim5 : read.index(adapter[:min_5match]) - trim3]
                quality = loaded_fq[read_no + 2]
                trimmed_quality = quality[
                    trim5 : read.index(adapter[:min_5match]) - trim3
                ]
                trimmed_reads.append(
                        (loaded_fq[read_no - 1], trimmed_read, loaded_fq[read_no + 1], trimmed_quality)
                    )
            else:
                too_short_count += 1
        else:
            adapter_missing_count += 1
    reject_low_quality(trimmed_reads, error)
    print("Reads with missing adapter: {}".format(adapter_missing_count))
    print("Reads too short after trimming: {}".format(too_short_count))
    print(
        "Retained {} of {} reads ({}%)\n\n".format(
            len(trimmed_reads),
            len(loaded_fq[1::4]),
            len(trimmed_reads) * 100 / len(loaded_fq[1::4]),
        )
    )
    return trimmed_reads

def reject_low_quality(trimmed_reads, max_error):
    """Reject trimmed reads based on mean quality score."""
    remaining_reads=[]
    rejected = 0
    for i in trimmed_reads:
        if mean_error(i[3]) < max_error:
            remaining_reads.append(i)
        else:
            rejected += 1
    print("\nReads with a mean error above {}: {}.".format(max_error, rejected))
    return remaining_reads


def phred33_to_error(qual):
    """Convert a phred33 score to error."""
    return 10 ** (-(ord(qual) - 33) / 10.0)


def mean_error(quality):
    """Calculate the mean quality score for a read."""
    return sum([phred33_to_error(q) for q in quality]) / len(quality)


def write_fastx(trimmed_reads, out_fx, fa, gz):
    """Write the trimmed reads to a fasta file."""
    if not gz:
        with open(out_fx, "w") as f:
            write_reads(f, trimmed_reads, fa)
    else:
        with gzip.open(out_fx, "wt") as f:
            write_reads(f, trimmed_reads, fa)


def write_reads(f, trimmed_reads, fa=True):
    for i in trimmed_reads:
        if fa:
            f.write(">{}\n".format(i[0]))
            f.write("{}\n".format(i[1]))
        else:
            f.write("{}\n".format(i[0]))
            f.write("{}\n".format(i[1]))
            f.write("{}\n".format(i[2]))
            f.write("{}\n".format(i[3]))


def flex_adapter_trim(
    fq, adapter, out_fa, min_len=18, trim5=0, trim3=0, min_5match=8, error=0.1, fa=True, gz=False
):
    """Trim the adapters from the reads, 5' and 3' clip if needed, and write to a fasta or fastq file."""
    loaded_fq = load_fastq(fq)
    trimmed_reads = trim_adapter(
        loaded_fq, adapter, min_len, trim5, trim3, min_5match, error
    )

    write_fastx(trimmed_reads, out_fa, fa, gz)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Trim read adapters (and 5' and 3' bases if needed) from fastq file and write to a fasta file."
    )
    parser.add_argument("fastq", help="The fastq file to trim.")
    parser.add_argument("adapter", help="The adapter sequence to trim.")
    parser.add_argument(
        "output", help="The output fastx file. Must be a fasta or fastq file."
    )
    parser.add_argument(
        "-m",
        "--min_len",
        help="The minimum length of the trimmed read.",
        default=18,
        type=int,
    )
    parser.add_argument(
        "-t5",
        "--trim5",
        help="The number of bases to trim from the 5' end after adapter removal.",
        default=0,
        type=int,
    )
    parser.add_argument(
        "-t3",
        "--trim3",
        help="The number of bases to trim from the 3' end after adapter removal.",
        default=0,
        type=int,
    )
    parser.add_argument(
        "-m5",
        "--min_5match",
        help="The minimum number of bases to match at the 5' end of the adapter.",
        default=8,
        type=int,
    )
    parser.add_argument(
        "-gz", "--gzip", help="Compress the output file with gzip.", action="store_true"
    )
    parser.add_argument("-e", "--error", help="The maximum error rate.", default=0.1, type=float)
    args = parser.parse_args()
    if args.output.endswith(".fa") or args.output.endswith(".fasta"):
        fa = True
    elif args.output.endswith(".fq") or args.output.endswith(".fastq"):
        fa = False
    else:
        print("Output file must be a fasta or fastq file.")
        exit()
    flex_adapter_trim(
        args.fastq,
        args.adapter,
        args.output,
        args.min_len,
        args.trim5,
        args.trim3,
        args.min_5match,
        args.error,
        fa,
        args.gzip,
    )