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
