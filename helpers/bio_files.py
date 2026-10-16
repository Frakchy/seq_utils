from pathlib import Path


def read_top_descriptions(path: str) -> list[str]:
    """
    Extracts the Description of the first hit for every query
    from a BLAST text report.

    Arguments:
    path: str, path to the BLAST results txt file.

    Returns list[str], one description per query that has hits.
    """
    descriptions = []
    in_hits_section = False
    description_width = 0
    with Path(path).open(encoding="utf-8") as file:
        for line in file:
            if line.startswith("Sequences producing significant alignments"):
                in_hits_section = True
                description_width = 0
            elif in_hits_section and not description_width:
                if line.startswith("Description"):
                    description_width = line.index("Name")
            elif in_hits_section and line.strip():
                descriptions.append(line[:description_width].strip())
                in_hits_section = False
    return descriptions


def write_lines(path: str, lines: list[str]) -> None:
    """
    Writes lines to a file.

    Arguments:
    path: str, path to the output file.
    lines: list[str].
    """
    with Path(path).open("w", encoding="utf-8") as file:
        for line in lines:
            file.write(f"{line}\n")
