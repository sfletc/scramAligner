## The scramAligner pipeline

The scramAligner pipeline is designer to align small RNA reads to references sequences with no mis-matches.  There is no requirement for reference index building.  The pipeline is ideally suited to analysis of siRNAs generated from transgenes or viruses.  Alignment files are generated for discrete read lengths (e.g. 21, 22 and 24nt).  Along with the core aligner (scram2), a command line script for alignment plotting is provided.  

----

### Installation

TODO: Add installation instructions

### Usage

NAME:

scramAligner- a command line Golang-based small RNA exact matching aligner


DESCRIPTION:

The SCRAM Aligner is a Golang-based command line application that performs exact matching alignment of small RNA sequences against a reference FASTA file. It is a core component of the SCRAM pipeline for small RNA exact matching alignment.


SYNOPSIS:
```
scramAligner -r <reference> -f <fileSet> -l <length> -o <outFilePrefix> [flags]
```

FLAGS:

```
-r, --alignTo <path>
    Path to the FASTA reference file.

-f, --fastxSet <path>
    Comma-separated path to the read file or set of read file replicates (collapsed FASTA, FASTA, or FASTQ). GZIPped files must have .gz file extension.

-t, --readFileType <type>
    Read file type: cfa (collapsed FASTA), fa (FASTA), fq (FASTQ) (default: fq).

-l, --length <string>
    Comma-separated read (sRNA) lengths to align.

-o, --outFilePrefix <path>
    Path to the output file prefix (len.csv will be appended).

--noSplit
    Do not split alignment count for each read by the number of times it aligns.

--noNorm
    Do not normalize read counts by library size (i.e., reads per million reads).

--minLen <integer>
    Minimum read length to include for RPMR normalization (default: 18).

--maxLen <integer>
    Maximum read length to include for RPMR normalization (default: 32).

--minCount <float>
    Minimum read count for alignment and to include for RPMR normalization (default: 1.0).

--adapter <sequence>
    3' adapter sequence to trim - FASTA & FASTQ only (default: "nil").
```

EXAMPLE:
```
scramAligner -r reference.fa -f reads.fq.gz -l 21,22,24 -o out
```
-------

## scramAligner Plotting Tool

The scramAligner plotting tool provides you with an easy and effective way to visualize the abundance profiles from your alignment outputs. The tool takes in various arguments to tailor your visualization as per your needs. It also provides options to plot coverage, abundance, and error based on your inputs.

USAGE: 

The tool can be run from the command line as follows:

```
python plot_tool.py [align_prefix] [align_lens] [header] [options]
```

#### Required Arguments:

- `align_prefix`: This argument specifies the prefix of alignment files - this is typically the prefix of the output file from scramAligner.
- `align_lens`: This argument requires a comma-separated list of siRNA lengths to plot.
- `header`: The exact header of the reference sequence (excluding the '>') to plot.

#### Optional Arguments:


- `-s`, `--smoothing_window`: Set the smoothing window (default is 1).
- `-c`, `--coverage`: Include this flag to plot coverage.
- `-a`, `--abundance`: Include this flag to plot abundance.
- `-e`, `--error`: Include this flag to plot error.
- `-y`, `--ylim`: Set y-axis limit. This argument takes in 2 values in a tuple indicating the lower and upper limits of the y-axis.
- `-n`, `--no_save`: Include this flag to not save the plot.

Note: If you choose to plot abundance, coverage will automatically be plotted.

Example:

```
python plot_tool.py align_prefix align_lens "Plot Header" -s 3 -c -a -y (-100,100)
```

This command will plot both coverage and abundance with a smoothing window of size 3 and y-axis limit from -100 to 100. The resulting plot will be saved by default.

------

### License

BSD 3-Clause License

### Contact

For questions or comments, please contact: Stephen Fletcher (s.fletcher@uq.edu.au)