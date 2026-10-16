# seq_utils

A set of bioinformatics utilities for working with FASTQ, BLAST and GenBank files. Written in pure Python (standard library only), requires Python 3.10+.

## Structure

- `main.py` - entry point with the public functions.
- `helpers/fastq.py` - helper functions for `filter_fastq`.
- `helpers/bio_files.py` - helper functions for `parse_blast_output` and `select_genes_from_gbk_to_fasta`.

## Functions

- `filter_fastq` - filters FASTQ reads by GC content, length and mean quality (phred33). Reads are processed one at a time, the result is saved to the `filtered/` directory.
- `parse_blast_output` - extracts the best hit (first Description) for each query from a BLAST text report and saves the names, sorted alphabetically, one per line.
- `select_genes_from_gbk_to_fasta` - saves protein sequences of the genes located before and after the given genes of interest in a GenBank file to a FASTA file.
