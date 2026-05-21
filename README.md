#GRAPES2
NOTE: This project is still in early developing phase!

## Installation

```bash
python3 install.py
```

## Options

Single-exon CNV analysis is enabled by default for BED files with up to 10,000
targets after exon splitting. Larger assays skip single-exon CNV analysis by
default while keeping segmented CNV analysis enabled.

Use `--single_exon_cnv` to force single-exon CNV analysis, or
`--skip_single_exon_cnv` to force it off.

## License
TODO
