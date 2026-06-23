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

import os
import tempfile

import numpy as np
import pytest

from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.runPlateletSim import run_platelet_sim
from runscripts.manual.runSecondWave import run_second_wave
from wholecell.io.tablereader import TableReader


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
	These pin it — the SOCE-dependent plateau, and the *mechanism* behind the
	sustained plateau.

	Dolan (2014) drove the cell with a transient IP3 dose and modelled neither
	thromboxane nor granule secretion, so its cytosol recovers after the peak.
	Our model adds two autocrine positive-feedback loops Dolan omitted —
	TXA2 -> TP -> Gq and secreted ADP -> P2Y1 — which sustain Gq -> IP3 and
	hold the cell activated (the platelet "second wave"). The high sustained
	plateau is therefore a model *prediction* (extra biology), not a recovery
	defect: disable the loops and apply a transient (reversible ADP) stimulus —
	the Dolan-equivalent configuration — and IP3 recovers toward its 50 nM
	baseline, so the Ca-handling machinery (IP3R / SERCA / SOCE / PMCA / NCX)
	reproduces Dolan's recovery. See
	reports/lab-books/lab-book-2026-06-14-recovery-phase.md.
	"""

	def test_plateau_sustained_and_active(self, _recovery_traces):
		"""+Ca: the cytosolic plateau stays sustained and active (>200 nM, the
		qualitative Dolan "active" criterion) — guards against a recovery-phase
		collapse."""
		plus, _ = _recovery_traces
		assert float(plus[120:200].mean()) > 200.0

	def test_sustained_plateau_requires_extracellular_ca(self, _recovery_traces):
		"""Dolan's central result extended past the peak: the sustained plateau
		is SOCE-dependent — the +Ca plateau stays well above the EDTA plateau
		(model ~430 vs ~235 nM over 120-200 s; NCX now extrudes under EDTA)."""
		plus, edta = _recovery_traces
		assert float(plus[120:200].mean()) > float(edta[120:200].mean()) + 100.0

	def test_autocrine_loops_sustain_else_ip3_recovers(self):
		"""The sustained plateau is the autocrine "second wave", not a defect.

		Contrast the full v0.61 model (autocrine loops on, thrombin drive)
		against the Dolan-equivalent configuration (loops off, transient
		reversible ADP stimulus). Both activate strongly (IP3 rises well above
		the 50 nM rest); but the full model holds IP3 elevated through the
		recovery phase (the TXA2/ADP loops sustain Gq), while the
		Dolan-equivalent recovers IP3 toward baseline — demonstrating the
		recovery machinery is intact and the sustained plateau is the autocrine
		amplification Dolan did not model. (The cytosolic Ca2+ and DTS-store
		recovery lag IP3 by ~100 s — buffer-limited — so IP3 is the clean,
		fast readout here; see the lab book / figure for the slower pools.)
		"""
		full = _cal(_run(300), 'ip3_nM')
		dolan_eq = _cal(_run(300, thrombin_peak_nM=0.0, atp_ex_peak_uM=0.0,
			autocrine_adp_gain=0.0, cox1_factor=0.0), 'ip3_nM')
		# Both configurations genuinely activate.
		assert full.max() > 150.0
		assert dolan_eq.max() > 150.0
		# Full model: IP3 stays elevated through the recovery phase (second wave).
		assert full[-1] > 150.0
		# Dolan-equivalent: IP3 recovers toward its 50 nM baseline once the
		# autocrine amplifiers are removed.
		assert dolan_eq[-1] < 120.0
		assert dolan_eq[-1] < full[-1] - 80.0


@pytest.mark.slow
class TestMcuCouplingDirection:
	"""#76 Part 2 — end-to-end regression guard for the MCU-knockout result.
	These pin the ODE's *behaviour*; the listener factor tests only check the
	formula, and the byte-identical goldens / Dolan 5/5 are wild-type-only
	(factor ≡ 1 at WT, so they are blind to the coupling by construction)."""

	def test_ko_reduces_transient_and_gain_toggle_flips_sign(self):
		"""Coupling ON: MCU KO *reduces* the evoked cytosolic peak + AUC (the
		headline). Coupling OFF (gain=0, bare buffer loss): KO *raises* it — the
		discriminating signature a deleted/broken/inverted coupling would fail."""
		wt = _cal(_run(300, agonist_delay_s=60.0, mcu_vmax_scale=1.0), 'ca_cyt_nM')
		ko = _cal(_run(300, agonist_delay_s=60.0, mcu_vmax_scale=0.0), 'ca_cyt_nM')
		ko_off = _cal(_run(300, agonist_delay_s=60.0, mcu_vmax_scale=0.0,
			mito_coupling_gain=0.0), 'ca_cyt_nM')
		assert ko.max() < wt.max()              # KO peak below WT
		assert ko.max() < 0.92 * wt.max()       # ~18% lower (loose band)
		assert ko.sum() < wt.sum()              # AUC lower (1 s steps)
		assert ko_off.max() > wt.max()          # buffer-only raises KO (sign flip)

	def test_ko_preserves_resting_cytosolic_ca(self):
		"""The evoked-specific gate spares the resting state (measured ~5% drift)."""
		rest = dict(thrombin_peak_nM=0.0, adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
		wt = _cal(_run(60, mcu_vmax_scale=1.0, **rest), 'ca_cyt_nM')
		ko = _cal(_run(60, mcu_vmax_scale=0.0, **rest), 'ca_cyt_nM')
		assert abs(ko[-1] - wt[-1]) < 0.12 * wt[-1]
