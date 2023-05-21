## The scramAligner pipeline

The scramAligner pipeline is designer to align small RNA reads to references sequences with no mis-matches.  There is no requirement for reference index building.  The pipeline is ideally suited to analysis of siRNAs generated from transgenes or viruses.  Alignment files are generated for discrete read lengths (e.g. 21, 22 and 24nt).  Along with the core aligner (scram2), command line scripts for adapter trimming (adaptTrim) and alignment plotting are provided.  

----

### Installation

- TODO

### Usage

NAME
scram2- a command line Golang-based small RNA exact matching aligner

SYNOPSIS
scram2 align [OPTIONS]

DESCRIPTION
The SCRAM Aligner is a Golang-based command line application that performs exact matching alignment of small RNA sequences against a reference FASTA file. It is a core component of the SCRAM pipeline for small RNA exact matching alignment.

OPTIONS
-h, --help
Display the help message and exit.

```-r, --alignTo <path>
    Path to the FASTA reference file.

-1, --fastxSet1 <path>
    Comma-separated path to the read file set 1 (collapsed FASTA, FASTA, or FASTQ). GZIPped files must have .gz file extension.

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

EXAMPLES
scram2 align -r reference.fa -1 reads.fq -l 21,22,24 -o out

