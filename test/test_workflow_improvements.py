import json
import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

from modules import cluster, normalize, readcount, segment
from modules.utils import (
    atomic_output_path,
    stage_manifest_matches,
    update_stage_manifest,
    validate_delimited_file,
    validate_vcf_file,
)


class SampleStub(SimpleNamespace):
    def add(self, key, value):
        setattr(self, key, value)


def _write_targetdepth_outputs(command, sample_name="sample"):
    stage_dir = command[command.index("-o") + 1]
    output_name = command[command.index("-n") + 1]
    depth_header = f"chr\tstart\tend\texon\tgc\tmap\t{sample_name}\n"
    depth_row = "chr1\t100\t200\tGENE1\t50\t100\t42\n"
    with open(
        os.path.join(stage_dir, f"{output_name}.read.counts.bed"),
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(depth_header)
        handle.write(depth_row)

    if "-d" in command:
        with open(
            os.path.join(stage_dir, f"{output_name}.per.base.coverage.bed"),
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(depth_header)
            handle.write(depth_row)

    with open(
        os.path.join(stage_dir, "summary_metrics.log"),
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "SAMPLE\tTOTAL_READS\tREADS_ON_TARGET\tREADS_CHRX\t%ROI\t"
            "MEAN_COVERAGE\tMEAN_COUNTS\tMEAN_ISIZE\tSD_ISIZE\t"
            "MEAN_COVERAGE_X\tMEAN_COUNTS_X\n"
        )
        handle.write(
            f"{sample_name}.bam\t1000\t800\t50\t80.0\t42.0\t10\t200\t20\t21.0\t5\n"
        )


class ConditionalCoverageTestCase(unittest.TestCase):
    def _run_read_depth(self, single_exon_cnv):
        temporary_directory = tempfile.TemporaryDirectory()
        output_dir = temporary_directory.name
        sample = SampleStub(name="sample")
        analysis = {
            "output_name": "cohort",
            "output_dir": output_dir,
            "bam_dir": "/inputs/bams.list",
            "reference": "/inputs/reference.fa",
            "ready_bed": "/inputs/panel.bed",
            "threads": 2,
            "single_exon_cnv": single_exon_cnv,
        }

        def fake_run(command, **kwargs):
            _write_targetdepth_outputs(command)
            return SimpleNamespace(returncode=0)

        annotation_identity = lambda current, *_args: current
        with patch.object(readcount, "annotate_gc", side_effect=annotation_identity), patch.object(
            readcount, "annotate_mappability", side_effect=annotation_identity
        ), patch.object(readcount.subprocess, "run", side_effect=fake_run) as mocked_run:
            samples, result = readcount.extract_read_depth(
                [sample],
                analysis,
                {"targetdepth": "/tools/target_depth.py"},
                {},
            )

        return temporary_directory, sample, result, mocked_run.call_args.args[0]

    def test_exome_mode_skips_per_base_coverage(self):
        temporary_directory, sample, result, command = self._run_read_depth(False)
        try:
            self.assertNotIn("-d", command)
            self.assertIsNone(result["per_base_coverage"])
            self.assertEqual(sample.mean_coverage, 42.0)
            self.assertTrue(
                os.path.isfile(
                    os.path.join(temporary_directory.name, "cohort.read.counts.bed")
                )
            )
        finally:
            temporary_directory.cleanup()

    def test_single_exon_mode_keeps_per_base_coverage(self):
        temporary_directory, _sample, result, command = self._run_read_depth(True)
        try:
            self.assertIn("-d", command)
            self.assertTrue(os.path.isfile(result["per_base_coverage"]))
        finally:
            temporary_directory.cleanup()

    def test_normalization_skips_per_base_stage_for_exomes(self):
        with tempfile.TemporaryDirectory() as output_dir:
            bed = os.path.join(output_dir, "panel.bed")
            with open(bed, "w", encoding="utf-8") as handle:
                handle.write("chr1\t100\t200\tGENE1\n")
            raw_depth = os.path.join(output_dir, "cohort.read.counts.bed")
            raw_df = pd.DataFrame(
                {
                    "chr": ["chr1"],
                    "start": [100],
                    "end": [200],
                    "exon": ["GENE1"],
                    "gc": [50],
                    "map": [100],
                    "sample": [42],
                }
            )
            raw_df.to_csv(raw_depth, sep="\t", index=False)
            normalized_df = raw_df.copy()
            normalized_df["sample_normalized_final"] = 1.0
            sample = SampleStub(name="sample", ontarget_reads=800)
            analysis = {
                "bed": bed,
                "unified_raw_depth": raw_depth,
                "output_name": "cohort",
                "output_dir": output_dir,
                "use_baseline_db": False,
                "offtarget": False,
                "single_exon_cnv": False,
            }

            with patch.object(
                normalize, "normalize_exon_level", return_value=normalized_df
            ), patch.object(normalize, "normalize_per_base") as mocked_per_base:
                normalize.launch_normalization([sample], analysis, {})

            mocked_per_base.assert_not_called()
            self.assertIsNone(analysis["normalized_per_base"])


class ReferenceSelectionTestCase(unittest.TestCase):
    def test_raw_threshold_and_maximum_reference_count_are_enforced(self):
        names = ["sample", "ref1", "ref2", "ref3", "poor"]
        raw_values = {
            "sample": [1.0, 0.95, 0.90, 0.86, 0.80],
            "ref1": [0.95, 1.0, 0.92, 0.88, 0.75],
            "ref2": [0.90, 0.92, 1.0, 0.87, 0.74],
            "ref3": [0.86, 0.88, 0.87, 1.0, 0.73],
            "poor": [0.80, 0.75, 0.74, 0.73, 1.0],
        }
        raw = pd.DataFrame(raw_values, index=names)
        similarity = (raw + 1.0) / 2.0

        with tempfile.TemporaryDirectory() as output_dir:
            correlation_file = os.path.join(output_dir, "correlation.tsv")
            similarity.to_csv(correlation_file, sep="\t")
            sample = SampleStub(name="sample", analysis_json={})

            cluster.cluster_samples(
                correlation_file,
                [sample],
                {"use_baseline_db": False},
                min_correlation=0.85,
                min_refs=2,
                max_refs=2,
            )

            self.assertEqual([name for name, _ in sample.references], ["ref1", "ref2"])
            self.assertEqual(sample.analysis_json["reference_count"], 2)
            self.assertAlmostEqual(
                sample.analysis_json["mean_raw_reference_correlation"],
                0.925,
            )
            self.assertEqual(sample.analyzable, "True")

    def test_insufficient_high_quality_references_marks_sample_unanalyzable(self):
        similarity = pd.DataFrame(
            [[1.0, 0.96], [0.96, 1.0]],
            index=["sample", "ref1"],
            columns=["sample", "ref1"],
        )
        with tempfile.TemporaryDirectory() as output_dir:
            correlation_file = os.path.join(output_dir, "correlation.tsv")
            similarity.to_csv(correlation_file, sep="\t")
            sample = SampleStub(name="sample", analysis_json={})
            cluster.cluster_samples(
                correlation_file,
                [sample],
                {"use_baseline_db": False},
                min_correlation=0.85,
                min_refs=2,
                max_refs=3,
            )
            self.assertEqual(sample.analyzable, "False")


class AtomicAndSegmentationTestCase(unittest.TestCase):
    def test_atomic_output_preserves_previous_file_after_failure(self):
        with tempfile.TemporaryDirectory() as output_dir:
            output_file = os.path.join(output_dir, "stage.tsv")
            with open(output_file, "w", encoding="utf-8") as handle:
                handle.write("old\n")

            with self.assertRaises(RuntimeError):
                with atomic_output_path(output_file) as temporary_path:
                    with open(temporary_path, "w", encoding="utf-8") as handle:
                        handle.write("new\n")
                    raise RuntimeError("interrupted")

            with open(output_file, "r", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "old\n")

    def test_manifest_rejects_output_when_an_input_changes(self):
        with tempfile.TemporaryDirectory() as output_dir:
            input_file = os.path.join(output_dir, "input.tsv")
            output_file = os.path.join(output_dir, "output.tsv")
            with open(input_file, "w", encoding="utf-8") as handle:
                handle.write("input\n")
            with open(output_file, "w", encoding="utf-8") as handle:
                handle.write("output\n")

            metadata = {"parameter": 1}
            update_stage_manifest(
                output_dir,
                "example",
                [output_file],
                metadata=metadata,
                input_paths=[input_file],
            )
            self.assertTrue(
                stage_manifest_matches(
                    output_dir,
                    "example",
                    [output_file],
                    input_paths=[input_file],
                    metadata=metadata,
                )
            )

            with open(input_file, "a", encoding="utf-8") as handle:
                handle.write("changed\n")
            self.assertFalse(
                stage_manifest_matches(
                    output_dir,
                    "example",
                    [output_file],
                    input_paths=[input_file],
                    metadata=metadata,
                )
            )

    def test_invalid_truncated_table_is_rejected(self):
        with tempfile.TemporaryDirectory() as output_dir:
            output_file = os.path.join(output_dir, "truncated.tsv")
            with open(output_file, "w", encoding="utf-8") as handle:
                handle.write("chr\tstart\tend\n")
            self.assertFalse(
                validate_delimited_file(
                    output_file,
                    required_columns=["chr", "start", "end"],
                    min_data_rows=1,
                )
            )

    def test_vcf_requires_complete_header(self):
        with tempfile.TemporaryDirectory() as output_dir:
            valid_vcf = os.path.join(output_dir, "valid.vcf")
            with open(valid_vcf, "w", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.3\n")
                handle.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                )
            self.assertTrue(validate_vcf_file(valid_vcf))

            truncated_vcf = os.path.join(output_dir, "truncated.vcf")
            with open(truncated_vcf, "w", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.3\n")
            self.assertFalse(validate_vcf_file(truncated_vcf))

    def test_segmentation_has_no_per_target_console_output_and_is_resumable(self):
        class FakeHmm:
            def __init__(self, *_args, **_kwargs):
                pass

            def fit_dispersion(self, **_kwargs):
                return 0.1

            def forward(self):
                return None

            def decode(self):
                return np.array([2, 2]), [[0, 0, 0.9], [0, 0, 0.8]]

            def posterior_decoding(self):
                return np.array(
                    [
                        [0.03, 0.03, 0.90, 0.02, 0.02],
                        [0.05, 0.05, 0.80, 0.05, 0.05],
                    ]
                )

            def calculate_map(self):
                return None, [0.9, 0.8]

            def compute_log_likelihood(self):
                return [[-2.0, -1.0, -0.1], [-2.0, -1.0, -0.2]]

        with tempfile.TemporaryDirectory() as output_dir:
            ratio_file = os.path.join(output_dir, "sample.ratios.bed")
            with open(ratio_file, "w", encoding="utf-8") as handle:
                handle.write("chr\tstart\tend\texon\tgc\tmap\tsample_ratio\n")
                handle.write("chr1\t100\t200\tGENE1\t50\t100\t0.1\n")
                handle.write("chr1\t300\t400\tGENE2\t50\t100\t0.2\n")
            sample = SampleStub(
                name="sample",
                analyzable="True",
                ratio_file=ratio_file,
                sample_folder=output_dir,
            )

            with patch.object(segment, "CustomHMM", FakeHmm), patch.object(
                segment, "calculate_positional_mean_variance", return_value={}
            ), patch("builtins.print") as mocked_print:
                segment.custom_hmm_seg([sample], {"output_dir": output_dir})

            mocked_print.assert_not_called()
            self.assertTrue(os.path.isfile(sample.segment_file))
            self.assertTrue(os.path.isfile(sample.segment_extended_file))
            self.assertTrue(os.path.isfile(sample.segment_file_map))
            with open(sample.segment_file, "r", encoding="utf-8") as handle:
                segment_fields = handle.readline().rstrip("\n").split("\t")
            self.assertAlmostEqual(float(segment_fields[-1]), 0.85)
            with open(sample.segment_file_map, "r", encoding="utf-8") as handle:
                first_map_fields = handle.readline().rstrip("\n").split("\t")
            self.assertEqual(first_map_fields[-2], "2")
            self.assertAlmostEqual(float(first_map_fields[-1]), 0.90)

            class FailingHmm:
                def __init__(self, *_args, **_kwargs):
                    raise AssertionError("validated outputs should have been reused")

            with patch.object(segment, "CustomHMM", FailingHmm), patch.object(
                segment, "calculate_positional_mean_variance", return_value={}
            ):
                segment.custom_hmm_seg([sample], {"output_dir": output_dir})

            manifest_path = os.path.join(output_dir, ".grapes2.stage_state.json")
            with open(manifest_path, "r", encoding="utf-8") as handle:
                manifest = json.load(handle)
            self.assertIn("segmentation:sample", manifest["stages"])


if __name__ == "__main__":
    unittest.main()
