"""Quantitative validation targets for the v0.6 / v0.61 modules.

The existing module tests check *direction* (knockout = 0, baseline > 0; aspirin
abolishes TXA2; the closed loop sustains Ca2+). These add *magnitude bands* — the
modules must reproduce the literature-anchored size of each effect, not just its
sign — plus the Dolan *recovery* phase, which the 30-s peak/SOCE "5/5" criteria
do not score.

Bands are anchored to the model's measured values (2026-06-14) with generous
margins, and to the qualitative literature claims: Mundell 2006 ("rapid,
substantial" PKC-dependent P2Y1 desensitisation), Purvis 2008 (the PLCb brake
lowers IP3), and the autocrine second wave. Windows are the shortest that expose
each signal — these are slow (sim-running) tests.
"""

import json
import os
import tempfile

import numpy as np
import pytest

from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.runPlateletSim import run_platelet_sim
from runscripts.manual.runSecondWave import run_second_wave
from wholecell.io.tablereader import TableReader


_REPO_ROOT = os.path.abspath(
	os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
_DOLAN_REF = os.path.join(_REPO_ROOT, 'reports', 'data',
	'dolan-2014-fig4-reference.json')


def _run(length, **config_kwargs):
	"""Run one platelet sim into a temp dir; return its simOut path.

	``config_kwargs`` are RunConfig fields (agonist peaks, perturbation scales).
	"""
	config_kwargs.setdefault('ca_ex_mM', 1.2)
	td = tempfile.mkdtemp()
	paths = run_platelet_sim(td, length_sec=length, seed=0,
		log_to_shell=False, run_config=RunConfig(**config_kwargs))
	return paths['sim_out_dir']


def _cal(sim_out, col):
	return TableReader(os.path.join(sim_out, 'CalciumTrace')
		).readColumn(col).flatten()


def _sec(sim_out, col):
	return TableReader(os.path.join(sim_out, 'SecretionTrace')
		).readColumn(col).flatten()


@pytest.mark.slow
class TestPKCFeedbackMagnitude:
	"""v0.6 brakes reach the size the literature reports, not just nonzero."""

	def test_p2y1_desensitisation_is_substantial(self):
		"""ADP-only: PKC desensitises a large fraction of active P2Y1.

		Mundell 2006 / Nicholas 2023 report rapid, *substantial* PKC-dependent
		P2Y1 desensitisation in human platelets. ADP-only isolates the P2Y1 arm
		(PARs off). The desensitised fraction must reach >= 0.5 (model ~0.69 at
		120 s, ~0.8 by 200 s) — more than the >0 the perturbation smoke test
		asserts — and stay zero at rest.
		"""
		frac = _cal(_run(120, thrombin_peak_nM=0.0, atp_ex_peak_uM=0.0),
			'p2y1_desensitised_frac')
		assert frac[0] == 0.0
		assert frac.max() >= 0.5

	def test_plcb_brake_holds_ip3_below_unbraked(self):
		"""Standard transient: the PLCb brake holds IP3 well below the
		no-brake (knockout) case — the Purvis 2008 result. Baseline IP3 is
		>= 15% below the unbraked plateau (model ~25%); resting IP3 unchanged.
		"""
		base = _cal(_run(150), 'ip3_nM')
		ko = _cal(_run(150, k_plcb_phos_scale=0.0), 'ip3_nM')  # PLCβ brake off
		assert abs(base[0] - ko[0]) < 1.0          # same resting IP3
		assert ko[-1] > base[-1]                    # brake lowers IP3
		assert base[-1] <= 0.85 * ko[-1]           # by >= 15%


@pytest.mark.slow
class TestSecretionKinetics:
	"""v0.61 secretion: dense-granule release leads alpha, and is substantial."""

	def test_dense_leads_alpha_and_release_is_substantial(self):
		"""Standard activation. Serotonin (dense) releases faster and further
		than fibrinogen (alpha) — dense granules lead (k_sec dense > alpha).
		Release is substantial and P-selectin is surface-exposed (the canonical
		flow-cytometry activation marker). Serotonin is the clean dense readout:
		it has no extracellular clearance, whereas secreted ADP is hydrolysed by
		ecto-NTPDase, which depresses its *released* fraction.
		"""
		sim = _run(120)
		sht = _sec(sim, 'serotonin_released_frac')  # dense marker
		fga = _sec(sim, 'fibrinogen_released_frac')  # alpha
		psel = _sec(sim, 'pselectin_surface_frac')   # alpha surface marker
		assert sht[0] == 0.0 and fga[0] == 0.0      # zero at rest
		assert sht[-1] > fga[-1]                     # dense leads alpha
		assert sht[-1] > 0.5                         # substantial dense release
		assert psel[-1] > 0.3                        # P-selectin exposed
		assert np.all(np.diff(sht) >= -1e-9)        # cargo only leaves


@pytest.mark.slow
class TestSecondWaveMagnitude:
	"""v0.61 capstone: the autocrine second wave reaches the reported size,
	and the decomposition holds (autocrine ADP carries Ca2+; TXA2 adds IP3)."""

	def test_sustained_gap_and_decomposition(self, tmp_path):
		out = str(tmp_path / 'sw')
		os.makedirs(out)
		_results, scalars = run_second_wave(out, adp_uM=0.5, length=300,
			log_to_shell=False)

		# The sustained second wave: closed loop vs open loop at 300 s
		# (model ~+93 nM).
		gap = scalars['full']['end_cyt_nM'] - scalars['v06']['end_cyt_nM']
		assert 40.0 <= gap <= 160.0
		# Autocrine ADP carries the sustained Ca2+ (aspirin ~ full on Ca2+) ...
		assert abs(scalars['aspirin']['end_cyt_nM']
			- scalars['full']['end_cyt_nM']) < 5.0
		# ... while TXA2 adds Gq drive seen in IP3, not Ca2+.
		assert (scalars['full']['ip3_end_nM']
			> scalars['aspirin']['ip3_end_nM'] + 20.0)


@pytest.fixture(scope='module')
def _recovery_traces():
	"""+Ca and EDTA cytosolic traces over a 200-s recovery window (shared)."""
	plus = _cal(_run(200), 'ca_cyt_nM')
	edta = _cal(_run(200, ca_ex_mM=0.0), 'ca_cyt_nM')
	return plus, edta


@pytest.mark.slow
class TestDolanRecoveryPhase:
	"""The 30-s peak/SOCE "5/5" criteria don't score the recovery phase.
	These pin it — including the SOCE-dependent plateau (passes) and a KNOWN
	GAP vs Dolan's quantitative plateau band (xfail)."""

	def test_plateau_sustained_and_active(self, _recovery_traces):
		"""+Ca: the cytosolic plateau stays sustained and active (>200 nM, the
		qualitative Dolan "active" criterion) — guards against a recovery-phase
		collapse."""
		plus, _ = _recovery_traces
		assert float(plus[120:200].mean()) > 200.0

	def test_sustained_plateau_requires_extracellular_ca(self, _recovery_traces):
		"""Dolan's central result extended past the peak: the sustained plateau
		is SOCE-dependent — the +Ca plateau stays well above the EDTA plateau
		(model ~430 vs ~250 nM over 120-200 s)."""
		plus, edta = _recovery_traces
		assert float(plus[120:200].mean()) > float(edta[120:200].mean()) + 100.0

	@pytest.mark.xfail(strict=True, reason=(
		"Known gap: the sustained +Ca plateau (~430 nM) sits above Dolan's "
		"200-275 nM plateau band under the saturating default agonist; the "
		"recovery phase is not yet calibrated (see the version-comparison doc "
		"validation section). Strict xfail flags it if a future fix closes the "
		"gap."))
	def test_plateau_matches_dolan_band(self, _recovery_traces):
		plus, _ = _recovery_traces
		with open(_DOLAN_REF) as f:
			ref = json.load(f)
		lo = ref['with_extracellular_ca']['plateau_cyt_nM']['range_lo']
		hi = ref['with_extracellular_ca']['plateau_cyt_nM']['range_hi']
		assert lo <= float(plus[120:200].mean()) <= hi
