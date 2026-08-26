# GRAPES2

GRAPES2 is a read-depth based CNV (Copy Number Variant) and SV (Structural
Variant / breakpoint) caller designed for targeted NGS sequencing — gene
panels, exomes, and similar amplicon/hybrid-capture designs. It takes a set
of aligned BAM files and a BED file of target regions, and outputs per-sample
CNV calls (BED/VCF) together with genome-wide log2-ratio plots.

> **Status:** early development. Interfaces, output formats and defaults may
> still change between releases.

## How it works

GRAPES2 runs a per-cohort pipeline over all BAM files given as input:

1. **Read depth extraction** — on-target (and optionally off-target) read
   counts are extracted per region using the bundled `TargetDepth` /
   `grapes_sv` C++ tools.
2. **Normalization** — read counts are normalized for GC content and
   mappability (using the bundled mappability/blacklist annotation tracks).
3. **Sample clustering** — samples are clustered/correlated to build an
   appropriate in-cohort or database-backed reference/baseline for each
   sample.
4. **Ratio calculation** — per-target log2 ratios are computed against the
   reference baseline.
5. **Segmentation** — a custom Hidden Markov Model (`modules/hmmnb.py`)
   segments the log2-ratio signal into copy-number states per sample.
6. **CNV calling** — segmented (multi-exon) and single-exon candidate calls
   are generated, scored and filtered (z-score, coefficient of variation,
   Gaussian/KL-divergence based likelihood).
7. **SV breakpoint calling** *(optional, `--breakpoint`)* — split-read/
   discordant-pair based structural variant detection.
8. **Export** — per-sample CNV/SV calls are merged into a VCF, a JSON
   summary, genome-wide CNV plots, and a cohort-wide "all calls" BED file.

## Installation

### Requirements

- Linux, Python 3.10
- `bedtools`, `samtools` available on `PATH`
- A C++ build toolchain (`build-essential`, `zlib1g-dev`, BLAS/LAPACK) if
  building the bundled C++ tools (`bin/grapes_sv`, `bin/htslib`,
  `bin/SeqLib`, `bin/TargetDepth`) from source

### Using Docker (recommended)

```bash
docker build -t grapes2 -f docker/dockerfile .
docker run --rm -v /path/to/data:/data grapes2 GRAPES2 --help
```

### Manual installation

```bash
git clone https://github.com/bdolmo/GRAPES2.git
cd GRAPES2
pip3 install -r requirements.txt
python3 install.py     # downloads mappability annotation tracks
```

`install.py` downloads the mappability tracks used for normalization into
`annotations/mappability/`. Blacklist regions and chromosome definitions for
`hg19`/`hg38` are already bundled under `annotations/`.

## Usage

```bash
python3 grapes2.py \
  --bam_dir /path/to/bams \
  --bed targets.bed \
  --output_dir /path/to/output \
  -f /path/to/reference.fasta \
  -g hg19
```

`--bam_dir` can be either a directory containing `.bam` files, or a plain
text file listing one BAM path per line (see `samples.list.txt` for an
example). All BAM files in a single run are treated as one cohort and used
to build each other's reference baseline unless `--use_baseline_db` /
`--baseline_db` is given.

### Key arguments

| Argument | Description |
|---|---|
| `--bam_dir` | Input directory (or list file) of BAM files **(required)** |
| `--bed` | BED file of target regions **(required)** |
| `--output_dir` | Output directory **(required)** |
| `-f`, `--fasta` | Reference genome FASTA **(required)** |
| `-g`, `--genome_version` | `hg19` or `hg38` **(required)** |
| `-t`, `--threads` | Number of CPU threads (default: `4`) |
| `--offtarget` | Also extract and analyze off-target coverage |
| `--breakpoint` | Perform SV breakpoint (split-read) analysis |
| `--single_exon_cnv` / `--skip_single_exon_cnv` | Force single-exon CNV reanalysis on/off. By default it's enabled for BED files with up to `--single_exon_cnv_target_limit` targets (10,000) after splitting, and skipped above that (e.g. exome BEDs) — see [Options](#options) below |
| `--use_baseline_db` / `--baseline_db` | Use a persistent SQLite baseline database instead of an in-cohort reference |
| `--upper_del_cutoff` | log2 ratio cutoff to call a deletion (default: `-0.6`) |
| `--lower_dup_cutoff` | log2 ratio cutoff to call a duplication (default: `0.4`) |
| `--min_zscore` | Minimum z-score to keep a call (default: `2.58`) |
| `--min_size` | Minimum reported CNV/SV size, in target count (default: `10`) |
| `--min_gc` / `--max_gc` | GC-content bounds for filtering targets (default: `20` / `80`) |
| `--min_mappability` | Minimum mappability for filtering targets (default: `30`) |
| `--plot_gene GENE1,GENE2` | Plot per-exon log2 ratios for the given gene(s) |
| `--force` | Force re-computation, ignoring cached intermediate files |

Run `python3 grapes2.py --help` for the full, up-to-date list of options.

## Output

For each sample `<sample>`, GRAPES2 writes to `<output_dir>/<sample>/`:

- `<sample>.GRAPES2.cnv.bed` — filtered/scored CNV calls
- `<sample>.GRAPES2.breakpoints.bed` — SV breakpoint calls (if `--breakpoint`)
- `<sample>.GRAPES2.bed` — merged CNV + SV calls
- `<sample>.GRAPES2.vcf` — final calls in VCF format
- `<sample>.GRAPES2.json` — structured per-sample analysis summary
- Genome-wide log2-ratio CNV plot (PNG)

At the cohort level, `<output_dir>/<output_dir_name>.all.calls.bed`
aggregates all samples' calls into a single file.

## Repository layout

```
grapes2.py          Command-line entry point / pipeline orchestration
modules/             Python pipeline modules (normalization, segmentation,
                      calling, plotting, VCF export, baseline DB, ...)
bin/                 Bundled C++ tools (TargetDepth, grapes_sv, htslib, SeqLib, samtools)
annotations/         Mappability, blacklist and chromosome reference tracks
ml_models/           Trained random-forest model used for variant scoring
docker/              Dockerfile for a containerized build
install.py           Post-clone setup: downloads annotation tracks
test/                Unit tests
```

## Options

Single-exon CNV analysis is enabled by default for BED files with up to 10,000
targets after exon splitting. Larger assays skip single-exon CNV analysis by
default while keeping segmented CNV analysis enabled.

Use `--single_exon_cnv` to force single-exon CNV analysis, or
`--skip_single_exon_cnv` to force it off.

## License

TODO
