from pathlib import Path
from typing import Iterator, TextIO


def normalize_bounds(bounds: float | tuple[float, float]) -> tuple[float, float]:
    """
    Converts bounds argument to a (lower, upper) tuple.

    Arguments:
    bounds: int / float / tuple of two numbers.

    Returns tuple (lower, upper).
    """
    if isinstance(bounds, (int, float)):
        return 0, bounds
    lower, upper = bounds
    return lower, upper


def make_output_path(output_fastq: str, output_dir: str = "filtered") -> Path:
    """
    Creates the output directory (if missing) and builds the output path.

    Arguments:
    output_fastq: str, name of the output file.
    output_dir: str, output directory.

    Returns Path to the output file.
    """
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    return directory / output_fastq


def read_fastq(path: str) -> Iterator[tuple[str, str, str, str]]:
    """
    Reads a FASTQ file, one record at a time.

    Arguments:
    path: str, path to the input FASTQ file.

    Yields tuples (header, sequence, separator, quality).
    """
    with Path(path).open() as file:
        while True:
            header = file.readline().rstrip("\n")
            if not header:
                break
            sequence = file.readline().rstrip("\n")
            separator = file.readline().rstrip("\n")
            quality = file.readline().rstrip("\n")
            yield header, sequence, separator, quality


def check_length(sequence: str, length_bounds: tuple[float, float]) -> bool:
    """
    Checks that sequence length lies within bounds.

    Arguments:
    sequence: str.
    length_bounds: tuple (lower, upper).

    Returns bool.
    """
    lower, upper = length_bounds
    return lower <= len(sequence) <= upper


def calculate_gc_content(sequence: str) -> float:
    """
    Computes GC content of a sequence in percent.

    Arguments:
    sequence: str.

    Returns float.
    """
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) * 100


def check_gc(sequence: str, gc_bounds: tuple[float, float]) -> bool:
    """
    Checks that GC content lies within bounds.

    Arguments:
    sequence: str, nucleotide sequence.
    gc_bounds: tuple (lower, upper) in percent.

    Returns bool.
    """
    lower, upper = gc_bounds
    return lower <= calculate_gc_content(sequence) <= upper


def calculate_mean_quality(quality: str) -> float:
    """
    Computes mean phred33 quality of a read.

    Arguments:
    quality: str, quality string of a read.

    Returns float.
    """
    if not quality:
        return 0.0
    total = sum(ord(symbol) - 33 for symbol in quality)
    return total / len(quality)


def check_quality(quality: str, quality_threshold: float) -> bool:
    """
    Checks that mean read quality is not below the threshold.

    Arguments:
    quality: str, quality string of a read (phred33).
    quality_threshold: int / float.

    Returns bool.
    """
    return calculate_mean_quality(quality) >= quality_threshold


def write_record(file: TextIO, record: tuple[str, str, str, str]) -> None:
    """
    Writes one FASTQ record (4 lines) to a file.

    Arguments:
    file: TextIO, file opened for writing.
    record: tuple (header, sequence, separator, quality).
    """
    file.write("\n".join(record) + "\n")
