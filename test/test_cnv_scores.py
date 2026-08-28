import math
import os
import tempfile
import unittest
from types import SimpleNamespace

from modules.call import (
    compute_sample_quality,
    export_cnv_calls_to_bed,
    score_multiple_exon,
    score_single_exon,
)
from modules.vcf import export_vcf_calls_to_bed


class CnvQualityTestCase(unittest.TestCase):
    def test_single_exon_returns_named_bounded_components(self):
        components = score_single_exon(
            "DEL",
            posterior_prob=0.9,
            sampleStd=0.12,
            pctCalls=1.0,
            corr=0.92,
            log2ratio=-1.0,
            control_cv=0.08,
            cn=1,
            nRois=1,
            kl_divergence=2.0,
            per_base_log2_snr=0.8,
            per_base_state_prob=0.9,
            return_components=True,
        )

        expected = {
            "CNV_QUALITY",
            "HMM_POSTERIOR",
            "SIGNAL_FIT",
            "DISPERSION_SCORE",
            "SAMPLE_QUALITY",
            "ROI_SUPPORT",
            "PERBASE_SUPPORT",
        }
        self.assertEqual(set(components), expected)
        for value in components.values():
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)

    def test_multi_exon_quality_uses_hmm_posterior(self):
        common = dict(
            svtype="DUP",
            sampleStd=0.12,
            pctCalls=1.0,
            corr=0.92,
            log2ratio=math.log2(1.5),
            event_std=0.08,
            cn=3,
            nRois=3,
        )
        low = score_multiple_exon(posterior_prob=0.1, **common)
        high = score_multiple_exon(posterior_prob=0.9, **common)
        self.assertGreater(high, low)
        self.assertAlmostEqual(high - low, 0.2)

    def test_roi_support_is_monotonic_without_eight_target_jump(self):
        common = dict(
            svtype="DUP",
            posterior_prob=0.9,
            sampleStd=0.12,
            pctCalls=1.0,
            corr=0.92,
            log2ratio=math.log2(1.5),
            event_std=0.08,
            cn=3,
            return_components=True,
        )
        two = score_multiple_exon(nRois=2, **common)
        eight = score_multiple_exon(nRois=8, **common)
        twenty = score_multiple_exon(nRois=20, **common)

        self.assertLess(two["ROI_SUPPORT"], eight["ROI_SUPPORT"])
        self.assertLess(eight["ROI_SUPPORT"], twenty["ROI_SUPPORT"])
        self.assertLess(eight["ROI_SUPPORT"], 1.0)
        self.assertLess(eight["CNV_QUALITY"], 1.0)

    def test_deletion_and_duplication_use_positive_symmetric_dispersion(self):
        common = dict(
            posterior_prob=0.9,
            sampleStd=0.12,
            pctCalls=1.0,
            corr=0.92,
            event_std=0.10,
            nRois=3,
        )
        deletion = score_multiple_exon(
            svtype="DEL",
            log2ratio=-1.0,
            cn=1,
            **common,
        )
        duplication = score_multiple_exon(
            svtype="DUP",
            log2ratio=math.log2(1.5),
            cn=3,
            **common,
        )
        self.assertAlmostEqual(deletion, duplication)

    def test_sample_quality_uses_raw_reference_correlation(self):
        low = compute_sample_quality(0.12, 1.0, 0.85)
        high = compute_sample_quality(0.12, 1.0, 0.95)
        self.assertGreater(high, low)

    def test_export_emits_quality_components_and_compatibility_alias(self):
        with tempfile.TemporaryDirectory() as output_dir:
            sample_dir = os.path.join(output_dir, "sample")
            os.mkdir(sample_dir)
            target_bed = os.path.join(output_dir, "targets.bed")
            with open(target_bed, "w", encoding="utf-8") as handle:
                handle.write("chr1\t100\t200\tGENE1\n")
                handle.write("chr1\t200\t300\tGENE2\n")
                handle.write("chr1\t300\t400\tGENE3\n")

            calls_bed = os.path.join(sample_dir, "sample.calls.bed")
            with open(calls_bed, "w", encoding="utf-8") as handle:
                handle.write(
                    "chr\tstart\tend\tregions\tgc\tmap\tz_score\tn_regions\t"
                    "log2_ratio\tcopy_number\thmm_posterior\tperbase_metrics\t"
                    "cnvtype\tcontrol_cv\tevent_std\n"
                )
                handle.write(
                    "chr1\t100\t400\tGENE1,GENE2,GENE3\t50\t100\t3.5\t3\t"
                    f"{math.log2(1.5):.6f}\t3\t0.9\t.\tDUP\t.\t0.08\n"
                )

            sample = SimpleNamespace(
                name="sample",
                analyzable="True",
                cnv_calls_bed=calls_bed,
                std_log2_ratio=0.12,
                log2_mad=0.12,
                mean_correlation=0.96,
                mean_raw_reference_correlation=0.92,
            )
            export_cnv_calls_to_bed(
                [sample],
                {"bed": target_bed, "output_dir": output_dir},
            )

            final_bed = os.path.join(
                sample_dir,
                "sample.GRAPES2.cnv.bed",
            )
            with open(final_bed, "r", encoding="utf-8") as handle:
                fields = handle.readline().rstrip("\n").split("\t")
            info = fields[3]
            self.assertIn("HMM_POSTERIOR=0.9", info)
            self.assertIn("QUALITY_MODEL=GRAPES2_HEURISTIC_V2", info)
            self.assertIn("DISPERSION_SCORE=", info)
            self.assertIn("CNV_QUALITY=", info)
            self.assertIn("CNV_SCORE=", info)

            info_values = {}
            for token in info.split(";"):
                if "=" in token:
                    key, value = token.split("=", 1)
                    info_values[key] = value
            self.assertEqual(info_values["CNV_QUALITY"], info_values["CNV_SCORE"])

    def test_passing_bed_excludes_explicitly_filtered_vcf_records(self):
        with tempfile.TemporaryDirectory() as output_dir:
            input_vcf = os.path.join(output_dir, "calls.vcf")
            output_bed = os.path.join(output_dir, "passing.bed")
            with open(input_vcf, "w", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.3\n")
                handle.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                )
                handle.write(
                    "chr1\t100\t.\tN\t<DUP>\t.\tPASS\tEND=200;CNV_QUALITY=0.9\tGT\t0/1\n"
                )
                handle.write(
                    "chr1\t300\t.\tN\t<DEL>\t.\tLow_RF_Score\tEND=400;CNV_QUALITY=0.8\tGT\t0/1\n"
                )

            export_vcf_calls_to_bed(
                input_vcf,
                output_bed,
                passing_only=True,
            )
            with open(output_bed, "r", encoding="utf-8") as handle:
                lines = handle.readlines()
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0].startswith("chr1\t100\t200\t"))


if __name__ == "__main__":
    unittest.main()
