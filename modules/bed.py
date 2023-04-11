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
