"""Byte-identical regression tests for the platelet simulation.

These tests freeze the *exact* numerical output of a small set of
deterministic scenarios. They exist to catch silent floating-point
drift introduced by refactors that aren't supposed to change behaviour
— most importantly the planned kinetics-as-data refactor (issue #32),
which has byte-identical Phase 3 output as an acceptance criterion.

The golden baselines are saved as NPZ files next to this test (under
`golden/`), and each scenario's CalciumTrace output is compared
column-by-column with `np.testing.assert_array_equal`. If the test
fails, the diff is inspectable: load the failing run's NPZ and the
golden NPZ in a notebook and plot side-by-side.

**Platform sensitivity.** Byte-identical equality is sensitive to
NumPy / SciPy version and to floating-point hardware. CI failures
caused by upstream library upgrades (not by repo changes) should be
investigated and, if benign, the golden NPZ regenerated with
`REGEN_GOLDEN=1 pytest models/platelet/tests/sim/test_byte_identical.py`
and the new file committed.
"""

from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np

from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.io.tablereader import TableReader


GOLDEN_DIR = os.path.join(os.path.dirname(__file__), 'golden')
REGEN = os.environ.get('REGEN_GOLDEN') == '1'

# Columns hashed for each scenario. Chosen to span the ODE state space:
# cytosolic Ca²⁺ (output), DTS Ca²⁺ (reservoir), IP3 (cascade output),
# SOCE flux (plasma-membrane influx).
TRACE_COLUMNS = ('ca_cyt_nM', 'ca_dts_uM', 'ip3_nM', 'soce_flux_nMs')


def _harvest_trace(sim_out_dir: str) -> dict[str, np.ndarray]:
	"""Read the CalciumTrace listener and return the byte-identical-able
	float64 columns as a dict keyed by column name."""
	reader = TableReader(os.path.join(sim_out_dir, 'CalciumTrace'))
	return {col: np.ascontiguousarray(
		reader.readColumn(col).flatten(), dtype=np.float64) for col in TRACE_COLUMNS}


def _run_scenario(scenario_name: str, **sim_kwargs) -> dict[str, np.ndarray]:
	"""Run one scenario and harvest its trace columns. The sim is run
	into a temp dir; only the trace returned to the caller is retained."""
	with tempfile.TemporaryDirectory() as td:
		paths = run_platelet_sim(td, log_to_shell=False, **sim_kwargs)
		return _harvest_trace(paths['sim_out_dir'])


def _golden_path(scenario_name: str) -> str:
	return os.path.join(GOLDEN_DIR, f'{scenario_name}.npz')


def _compare_or_regenerate(scenario_name: str,
		current: dict[str, np.ndarray]) -> None:
	"""Compare current trace to the golden NPZ, or regenerate if asked.

	Set REGEN_GOLDEN=1 in the environment to bootstrap or refresh the
	golden file. Without that, a missing golden file fails the test
	loudly — never silently captures a new baseline.
	"""
	path = _golden_path(scenario_name)
	if REGEN:
		os.makedirs(GOLDEN_DIR, exist_ok=True)
		np.savez(path, **current)
		raise unittest.SkipTest(
			f'REGEN_GOLDEN=1: wrote new baseline to {path}; '
			f'commit the file and re-run without the env var.')
	if not os.path.isfile(path):
		raise AssertionError(
			f'Golden baseline missing for scenario {scenario_name!r} at '
			f'{path}. Run with REGEN_GOLDEN=1 to capture the current '
			f'output as the new baseline.')
	golden = np.load(path)
	for col in TRACE_COLUMNS:
		np.testing.assert_array_equal(
			current[col], golden[col],
			err_msg=f'Byte-identical regression in scenario '
				f'{scenario_name!r}, column {col!r}. If this drift is '
				f'expected (e.g. you intentionally changed kinetics), '
				f'regenerate with REGEN_GOLDEN=1.')


class TestByteIdenticalRegression(unittest.TestCase):
	"""Locks in the bit-level numerical output of a small set of scenarios.

	Designed to bracket the kinetics-as-data refactor (#32) so that a
	parameter-externalisation step that is supposed to be behaviour-
	preserving can be verified byte-identical.
	"""

	def test_default_activation_30s(self):
		"""Default agonist stimulation (+Ca²⁺), 30 s. Exercises the full
		GPCR cascade → IP3 → IP3R → DTS-release → SOCE → PMCA loop."""
		trace = _run_scenario(
			'default_activation_30s',
			length_sec=30, seed=0, ca_ex_mM=1.2)
		_compare_or_regenerate('default_activation_30s', trace)

	def test_at_rest_30s(self):
		"""All agonist peaks zero, 30 s. Exercises the resting fixed
		point — receptors stay at REST, no IP3R firing, SOCE at basal."""
		trace = _run_scenario(
			'at_rest_30s',
			length_sec=30, seed=0, ca_ex_mM=1.2,
			thrombin_peak_nM=0.0, adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
		_compare_or_regenerate('at_rest_30s', trace)


if __name__ == '__main__':
	unittest.main()
