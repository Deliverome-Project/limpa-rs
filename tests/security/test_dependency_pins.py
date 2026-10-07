import copy
import importlib.util
import json
from pathlib import Path
import tomllib
import unittest

ROOT = Path(__file__).parents[2]
spec = importlib.util.spec_from_file_location('pins', ROOT / 'scripts/check-dependency-pins.py')
pins = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pins)


class PinPolicyTests(unittest.TestCase):
    def setUp(self):
        self.project = tomllib.loads((ROOT / 'pyproject.toml').read_text())
        self.python_lock = tomllib.loads((ROOT / 'uv.lock').read_text())
        self.r_lock = json.loads((ROOT / 'renv.lock').read_text())
        self.bootstrap = (ROOT / 'renv/activate.R').read_text()
        self.cargo = tomllib.loads((ROOT / 'Cargo.toml').read_text())
        self.cargo_lock = tomllib.loads((ROOT / 'Cargo.lock').read_text())

    def test_current_pins_pass(self):
        self.assertEqual(pins.check_python(self.project, self.python_lock), [])
        self.assertEqual(pins.check_r(self.r_lock, self.bootstrap), [])
        self.assertEqual(pins.check_rust(self.cargo, self.cargo_lock), [])

    def test_python_range_and_missing_artifact_hash_fail(self):
        self.project['project']['dependencies'][0] = 'numpy>=1.26'
        self.assertTrue(pins.check_python(self.project, self.python_lock))
        self.setUp()
        self.python_lock['package'][0]['wheels'][0].pop('hash')
        self.assertTrue(pins.check_python(self.project, self.python_lock))

    def test_python_pin_absent_from_lock_fails(self):
        self.project['build-system']['requires'][0] = 'setuptools==999.0.0'
        self.assertTrue(pins.check_python(self.project, self.python_lock))

    def test_r_missing_version_and_floating_revision_fail(self):
        for field, value in [('Version', 'latest'), ('RemoteSha', 'RELEASE_3_23')]:
            lock = copy.deepcopy(self.r_lock)
            lock['Packages']['limpa'][field] = value
            self.assertTrue(pins.check_r(lock, self.bootstrap))

    def test_r_transitive_dependency_and_bootstrap_drift_fail(self):
        del self.r_lock['Packages']['statmod']
        self.assertTrue(pins.check_r(self.r_lock, self.bootstrap))
        self.setUp()
        self.assertTrue(pins.check_r(self.r_lock, self.bootstrap + '\n# changed'))

    def test_rust_ranges_git_branches_and_missing_checksums_fail(self):
        for requirement in ['0.34', {'git': 'https://example.com/crate', 'branch': 'main'}]:
            manifest = copy.deepcopy(self.cargo)
            manifest['dependencies']['nalgebra'] = requirement
            self.assertTrue(pins.check_rust(manifest, self.cargo_lock))
        self.cargo_lock['package'][0].pop('checksum')
        self.assertTrue(pins.check_rust(self.cargo, self.cargo_lock))

    def test_install_bypasses_fail(self):
        commands = ['python -m pip install numpy', 'uv pip install --no-deps numpy',
                    'uv sync', 'uv run --with numpy analysis.py',
                    'install.packages("limpa")', 'BiocManager::install("limpa")',
                    'remotes::install_github("SmythLab/limpa")', 'renv::install("limpa")',
                    'cargo build', 'cargo test --release', 'python -m build']
        for command in commands:
            with self.subTest(command=command):
                self.assertTrue(pins.check_commands(command, 'workflow'))

    def test_locked_commands_pass(self):
        commands = ['uv sync --locked --group ci --no-install-project --no-build',
                    'uv pip install --no-deps dist/*.whl', 'renv::restore(prompt=FALSE)',
                    'cargo build --release --locked', 'python -m build --no-isolation']
        for command in commands:
            self.assertEqual(pins.check_commands(command, 'workflow'), [])

    def test_r_import_not_locked_fails(self):
        self.assertTrue(pins.check_r_imports('library(newpackage)', self.r_lock['Packages'], 'analysis.R'))
        self.assertTrue(pins.check_r_imports('newpackage::fit(x)', self.r_lock['Packages'], 'analysis.R'))
        self.assertEqual(pins.check_r_imports('library(limpa); stats::median(x)', self.r_lock['Packages'], 'analysis.R'), [])
