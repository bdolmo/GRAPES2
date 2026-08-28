import os
import tempfile
import unittest

from modules.utils import remove_tmp_files


class IntermediateCleanupTestCase(unittest.TestCase):
    def test_removes_intermediates_and_preserves_final_outputs(self):
        intermediate_names = [
            "sample.bam_counts.bed",
            "sample.bam_isizes.bed",
            "sample.bam_coverage.bed",
            "cohort.read.counts.bed",
            "cohort.per.base.coverage.bed",
            "cohort.normalized.depth.bed",
            "cohort.normalized.per.base.bed",
            "cohort.ratios.bed",
        ]
        final_names = [
            "cohort.all.calls.bed",
            "cohort.grapes.log",
            "sample.GRAPES2.bed",
            "summary_metrics.log",
        ]

        with tempfile.TemporaryDirectory() as output_dir:
            for filename in intermediate_names + final_names:
                with open(os.path.join(output_dir, filename), "wb") as handle:
                    handle.write(b"1234")

            sample_dir = os.path.join(output_dir, "sample")
            os.makedirs(sample_dir)
            sample_vcf = os.path.join(sample_dir, "sample.GRAPES2.vcf")
            sample_json = os.path.join(sample_dir, "sample.GRAPES2.json")
            for file_path in (sample_vcf, sample_json):
                with open(file_path, "wb") as handle:
                    handle.write(b"final")

            removed_files, reclaimed_bytes = remove_tmp_files(output_dir)

            self.assertEqual(
                {os.path.basename(path) for path in removed_files},
                set(intermediate_names),
            )
            self.assertEqual(reclaimed_bytes, 4 * len(intermediate_names))
            for filename in intermediate_names:
                self.assertFalse(os.path.exists(os.path.join(output_dir, filename)))
            for filename in final_names:
                self.assertTrue(os.path.isfile(os.path.join(output_dir, filename)))
            self.assertTrue(os.path.isfile(sample_vcf))
            self.assertTrue(os.path.isfile(sample_json))


if __name__ == "__main__":
    unittest.main()
