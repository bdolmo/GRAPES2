# GRAPES2

GRAPES2 is a read-depth based CNV (Copy Number Variant) and SV (Structural
Variant / breakpoint) caller designed for targeted NGS sequencing — gene
panels, exomes, and hybrid-capture designs. It takes a set
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
```

The image does not bundle the mappability annotation tracks (~460MB
combined for hg19+hg38) — mount a host directory with them instead. Populate
that directory once, outside Docker, with `python3 install.py` (see below),
then mount it read-only on every run:

```bash
docker run --rm \
  -v /path/to/mappability:/usr/src/app/GRAPES2/annotations/mappability:ro \
  -v /path/to/data:/data \
  grapes2 python grapes2.py --bam_dir /data/bams --bed /data/targets.bed \
    --output_dir /data/out -f /data/reference.fasta -g hg19
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
| `--min_reference_correlation` | Minimum raw Spearman correlation for reference selection (default: `0.85`) |
| `--min_reference_samples` / `--max_reference_samples` | Required and maximum selected references per sample (default: `3` / `10`) |
| `--use_baseline_db` / `--baseline_db` | Use a persistent SQLite baseline database instead of an in-cohort reference |
| `--upper_del_cutoff` | log2 ratio cutoff to call a deletion (default: `-0.6`) |
| `--lower_dup_cutoff` | log2 ratio cutoff to call a duplication (default: `0.4`) |
| `--min_zscore` | Minimum z-score to keep a call (default: `2.58`) |
| `--min_cnv_quality` | Minimum final uncalibrated CNV quality to report (default: `0.5`) |
| `--min_size` | Minimum reported CNV/SV size, in target count (default: `10`) |
| `--min_gc` / `--max_gc` | GC-content bounds for filtering targets (default: `20` / `80`) |
| `--min_mappability` | Minimum mappability for filtering targets (default: `30`) |
| `--plot_gene GENE1,GENE2` | Plot per-exon log2 ratios for the given gene(s) |
| `--force` | Force re-computation, ignoring cached intermediate files |
| `--keep_intermediate_files` | Keep large coverage/count/normalization intermediates; they are removed after a successful run by default |

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

### CNV quality fields

`CNV_QUALITY` is the canonical GRAPES2 call-quality value. It is an
uncalibrated heuristic score in the range 0–1, not a posterior probability or
a Phred score. `CNV_SCORE` is emitted with the same value as a deprecated alias
for compatibility with existing consumers. Calls are filtered using
`--min_cnv_quality` (default `0.5`), and `QUALITY_MODEL` identifies the scoring
schema (`GRAPES2_HEURISTIC_V2`).

The final BED and VCF INFO fields expose the components used to calculate it:
`HMM_POSTERIOR`, `SIGNAL_FIT`, `DISPERSION_SCORE`, `SAMPLE_QUALITY`,
`ROI_SUPPORT`, and (for single-exon calls) `PERBASE_SUPPORT`. Raw supporting
dispersion values are reported separately as `EVENT_STD` and `CONTROL_CV`.
When random-forest scoring is enabled, `RF_SCORE` remains a separate,
uncalibrated model score and is not combined numerically with `CNV_QUALITY`.
The VCF retains low-RF records with `FILTER=Low_RF_Score`; passing BED and
cohort outputs exclude records carrying an explicit VCF filter.

After those final outputs have been written successfully, GRAPES2 removes
reproducible intermediate coverage, count, normalization, ratio, and temporary
BED files by default. Pass `--keep_intermediate_files` when they are needed for
debugging or inspection. Intermediates are retained automatically when a run
fails before final output generation.

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
