from pathlib import Path
from typing import Iterable

from helpers.bio_files import (
    find_neighbor_indices,
    format_cds_header,
    normalize_genes,
    parse_gbk,
    read_top_descriptions,
    select_translated_cds,
    write_fasta_record,
    write_lines,
)
from helpers.fastq import (
    check_gc,
    check_length,
    check_quality,
    make_output_path,
    normalize_bounds,
    read_fastq,
    write_record,
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: float | tuple[float, float] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters FASTQ reads by GC content, length and mean quality.

    Arguments:
    input_fastq: str, path to the input FASTQ file.
    output_fastq: str, name of the output file.
    gc_bounds: number or tuple (lower, upper), GC content in percent.
    length_bounds: number or tuple (lower, upper), read length.
    quality_threshold: int / float, minimal mean phred33 quality.

    Returns None.
    """
    gc_bounds = normalize_bounds(gc_bounds)
    length_bounds = normalize_bounds(length_bounds)
    output_path = make_output_path(output_fastq)

    with output_path.open("w") as output_file:
        for record in read_fastq(input_fastq):
            _, sequence, _, quality = record
            if (
                check_length(sequence, length_bounds)
                and check_gc(sequence, gc_bounds)
                and check_quality(quality, quality_threshold)
            ):
                write_record(output_file, record)


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Takes the best hit for every query from a BLAST
    text report and saves them alphabetically sorted.

    Arguments:
    input_file: str, path to the BLAST results txt file.
    output_file: str, path to the output file.

    Returns None.
    """
    descriptions = read_top_descriptions(input_file)
    write_lines(output_file, sorted(descriptions))


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: str | Iterable[str],
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = "output.fasta",
) -> None:
    """
    Saves protein sequences of genes located near genes of interest to a FASTA file.

    Arguments:
    input_gbk: str, path to the input GBK file.
    genes: str or collection of str, gene names.
    n_before: int >= 1.
    n_after: int >= 1.
    output_fasta: str, path to the output FASTA file.

    Returns None.
    Raises ValueError if n_before or n_after is less than 1.
    """
    if n_before < 1 or n_after < 1:
        raise ValueError("n_before and n_after must be >= 1")

    genes = normalize_genes(genes)
    records = parse_gbk(input_gbk)
    all_cds = [select_translated_cds(record) for record in records]

    with Path(output_fasta).open("w", encoding="utf-8") as output_file:
        for cds_list in all_cds:
            indices = find_neighbor_indices(cds_list, genes, n_before, n_after)
            for index in indices:
                cds = cds_list[index]
                header = format_cds_header(cds)
                sequence = cds["qualifiers"]["translation"]
                write_fasta_record(output_file, header, sequence)
