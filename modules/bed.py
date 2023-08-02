from intervaltree import Interval, IntervalTree
from collections import defaultdict

class BedRecord:
    """
    Class to store BED records.
    """

    def __init__(self, chr: str, start: int, end: int):
        self.chr = chr
        self.start = start
        self.end = end
        if start > end:
            msg = f" ERROR: start {start} cannot be greater than end {end}"
            raise ValueError(msg)

    def __repr__(self):
        return f"{self.chr}\t{self.start}\t{self.end}"


def load_bed_file(bed_file):
    """
    """
    regions = defaultdict(IntervalTree)
    with open(bed_file) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            rec = BedRecord(tmp[0], int(tmp[1]), int(tmp[2]))
            regions[rec.chr][rec.start:rec.end] = rec
    return regions