import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


SCRIPTS = Path(__file__).parents[1] / "scripts"


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fit = load_script("fit_twas_weights")
assemble = load_script("assemble_twas_weight_list")
batches = load_script("twas_batches")


class TestTwasWeightWorker(unittest.TestCase):
    def test_adaptive_batch_ranges(self):
        self.assertEqual(batches.batch_ranges(0, 200), [])
        self.assertEqual(batches.batch_ranges(5, 200), [(1, 5)])
        self.assertEqual(batches.batch_ranges(400, 200), [(1, 200), (201, 400)])
        self.assertEqual(batches.batch_ranges(401, 200), [(1, 200), (201, 400), (401, 401)])

    def worker_fixture(self, root):
        args = SimpleNamespace(
            plink="plink",
            geno=root / "geno",
            covar=root / "covar.tsv",
            rscript="Rscript",
            fusion_script=root / "FUSION.compute_weights.R",
            gcta=root / "gcta",
            gemma=root / "gemma",
        )
        model = fit.Model(1, "1", 1_000_000, "gene1", ("2.5",))
        return args, model, root / "tmp", root / "weights"

    def test_missing_weight_classification(self):
        expected = [
            ("heritability 0.01; LRT P-value 0.5 : skipping gene", "insufficient heritability"),
            ("likely GCTA could not converge, skipping gene", "GCTA did not converge"),
            ("Only one SNP available, skipping this gene", "only one cis-SNP"),
        ]
        for message, diagnostic in expected:
            with self.subTest(message=message):
                self.assertEqual(fit.classify_missing_weight(message), ("skipped", diagnostic))
        self.assertEqual(
            fit.classify_missing_weight("ERROR: invalid covariate"),
            ("failed", "FUSION exited without producing a weight file"),
        )

    def test_read_exact_batch(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / "models.bed.gz"
            with gzip.open(bed, "wt") as handle:
                handle.write("#chr\tstart\tend\tphenotype_id\tS1\n")
                for index in range(1, 6):
                    handle.write(f"1\t{index}\t{index}\tgene{index}\t{index}.0\n")
            samples, models = fit.read_batch(bed, 2, 4)
            self.assertEqual(samples, ["S1"])
            self.assertEqual([model.bed_index for model in models], [2, 3, 4])
            self.assertEqual([model.phenotype_id for model in models], ["gene2", "gene3", "gene4"])

    def test_status_is_atomic(self):
        with tempfile.TemporaryDirectory() as directory:
            status = Path(directory) / "status.tsv"
            row = {field: "" for field in fit.STATUS_FIELDS}
            row.update(bed_index=1, phenotype_id="gene1", outcome="skipped")
            fit.write_status(status, [row])
            self.assertTrue(status.exists())
            self.assertFalse(Path(f"{status}.tmp").exists())

    def test_weight_assembly_uses_status_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            status = root / "status.tsv"
            rows = []
            for index, outcome, gene in [(3, "success", "g3"), (1, "success", "g1"), (2, "skipped", "g2")]:
                row = {field: "" for field in fit.STATUS_FIELDS}
                row.update(
                    bed_index=index,
                    phenotype_id=gene,
                    outcome=outcome,
                    weight_path=str(root / "expression" / f"{gene}.wgt.RDat"),
                )
                rows.append(row)
            fit.write_status(status, rows)
            self.assertEqual(
                assemble.successful_weights([status]),
                [root / "expression" / "g1.wgt.RDat", root / "expression" / "g3.wgt.RDat"],
            )

    def test_success_publishes_weight_and_cleans_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args, model, batch_temp, weights_dir = self.worker_fixture(root)

            def fake_run(command, log_handle):
                if "--save_hsq" in [str(item) for item in command]:
                    work_dir = Path(log_handle.name).parent
                    (work_dir / "weight.wgt.RDat").write_text("weight")
                    (work_dir / "weight.hsq").write_text("staged 0.2 0.05 0.001\n")

            with patch.object(fit, "run_logged", side_effect=fake_run):
                row = fit.fit_model(model, [("F1", "S1")], args, batch_temp, weights_dir)
            self.assertEqual(row["outcome"], "success")
            self.assertTrue((weights_dir / "gene1.wgt.RDat").exists())
            self.assertFalse(batch_temp.joinpath("000000001_gene1").exists())

    def test_skip_cleans_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args, model, batch_temp, weights_dir = self.worker_fixture(root)

            def fake_run(command, log_handle):
                if "--save_hsq" in [str(item) for item in command]:
                    log_handle.write("heritability 0.01; LRT P-value 0.5 : skipping gene\n")
                    log_handle.flush()

            with patch.object(fit, "run_logged", side_effect=fake_run):
                row = fit.fit_model(model, [("F1", "S1")], args, batch_temp, weights_dir)
            self.assertEqual(row["outcome"], "skipped")
            self.assertFalse(batch_temp.joinpath("000000001_gene1").exists())

    def test_empty_cis_window_skips_and_cleans_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args, model, batch_temp, weights_dir = self.worker_fixture(root)

            def fake_run(command, log_handle):
                log_handle.write("Error: All variants excluded.\n")
                log_handle.flush()
                raise fit.subprocess.CalledProcessError(12, command)

            with patch.object(fit, "run_logged", side_effect=fake_run):
                row = fit.fit_model(model, [("F1", "S1")], args, batch_temp, weights_dir)
            self.assertEqual(row["outcome"], "skipped")
            self.assertEqual(row["diagnostic"], "no variants in cis-window")
            self.assertFalse(batch_temp.joinpath("000000001_gene1").exists())

    def test_other_plink_exit_12_remains_fatal(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args, model, batch_temp, weights_dir = self.worker_fixture(root)

            def fake_run(command, log_handle):
                log_handle.write("Error: malformed genotype input.\n")
                log_handle.flush()
                raise fit.subprocess.CalledProcessError(12, command)

            with patch.object(fit, "run_logged", side_effect=fake_run):
                row = fit.fit_model(model, [("F1", "S1")], args, batch_temp, weights_dir)
            self.assertEqual(row["outcome"], "failed")
            self.assertTrue(batch_temp.joinpath("000000001_gene1", "model.log").exists())

    def test_unexpected_missing_weight_retains_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args, model, batch_temp, weights_dir = self.worker_fixture(root)

            def fake_run(command, log_handle):
                if "--save_hsq" in [str(item) for item in command]:
                    log_handle.write("unexpected model error\n")
                    log_handle.flush()

            with patch.object(fit, "run_logged", side_effect=fake_run):
                row = fit.fit_model(model, [("F1", "S1")], args, batch_temp, weights_dir)
            self.assertEqual(row["outcome"], "failed")
            self.assertTrue(batch_temp.joinpath("000000001_gene1", "model.log").exists())


if __name__ == "__main__":
    unittest.main()
