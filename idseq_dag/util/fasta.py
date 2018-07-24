#!/usr/bin/env python3
from typing import Iterator, List, Tuple, NamedTuple

class Read(NamedTuple):
    header: str
    sequence: str

def iterator(fasta_file: str) -> Iterator[Read]:
    """Iterate through fasta_file, yielding one Read tuple at a time."""
    # TODO: Support full fasta format, where sequences may be split over multiple lines.
    with open(fasta_file, 'r', encoding='utf-8') as f:
        while True:
            header = f.readline().rstrip()
            sequence = f.readline().rstrip()
            if not header or not sequence:
                break
            yield Read(header, sequence)

def synchronized_iterator(fasta_files: List[str]) -> Iterator[Tuple[Read, ...]]:
    """Iterate through one or more fasta files in lockstep, yielding tuples of
    matching reads.  When the given list fasta_files has length 1, yield
    1-element tuples.  This facilitates uniform processing of either
    unpaired or paired-end reads."""
    return zip(*map(iterator, fasta_files))
