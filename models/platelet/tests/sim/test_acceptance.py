"""Behavioural acceptance / regression suite for the platelet model.

This is the model's headline acceptance contract: a small set of *behavioural*
checks, each pinning one biological result the model must reproduce, written as
readable assertions with explicit, biologically-anchored tolerance bands.

It replaces two older, more brittle ideas:

  * the **byte-identical goldens** (the former ``test_byte_identical.py`` +
    ``golden/*.npz``) — bit-exact equality on four columns of two scenarios,
    sensitive to NumPy / SciPy versions and floating-point hardware, and blind
    to *meaning* (it could not say whether a drift mattered biologically); and
  * the **"Dolan 5/5" gate** (a single ``sum(criteria) == 5`` assertion) — an
    all-or-nothing meta-check that passed largely by construction (each new
    layer is normalised at rest to keep it green) and reported nothing about
    *which* behaviour moved.

Here each behaviour is its own test with its own band, so a regression names the
biology that broke. Bands are anchored to measured current values (seed 0,
2026-06-26) with margins chosen to catch a meaningful drift while tolerating the
run-to-run noise from ``RestingDecay`` (which draws from NumPy's global RNG); the
Ca²⁺/secretion/thromboxane/integrin trace columns are themselves deterministic at
seed 0.

Deeper, magnitude-level guards live in the subsystem suites
(``test_validation_targets.py``, ``test_thromboxane.py``, ``test_integrin.py``,
``test_inhibitory_axis.py``); structural invariants (mass > 0, SOCE ≥ 0, the
store never draining below cytosol) live in ``test_regression.py``. Together
these are the model's wider validation portfolio — see
``docs/validation-and-regressions.md``.

Every test here is ``slow`` (it runs the simulation). Run just this suite with::

    PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 \\
        python -m pytest models/platelet/tests/sim/test_acceptance.py -v
"""

import os
import tempfile

import pytest

from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.io.tablereader import TableReader


# All three receptors held at their REST level — an un-stimulated cell.
_AT_REST = dict(thrombin_peak_nM=0.0, adp_peak_uM=0.0, atp_ex_peak_uM=0.0)


def _run(length, **config_kwargs):
	"""Run one platelet sim into a temp dir; return its simOut path.

	``config_kwargs`` are :class:`RunConfig` fields (agonist peaks, knockout
	scales, …). ``ca_ex_mM`` defaults to the physiological 1.2 mM.
	"""
	config_kwargs.setdefault('ca_ex_mM', 1.2)
	td = tempfile.mkdtemp()
	paths = run_platelet_sim(td, length_sec=length, seed=0,
		log_to_shell=False, run_config=RunConfig(**config_kwargs))
	return paths['sim_out_dir']


def _col(sim_out, listener, column):
	"""Read one listener column from a finished run as a flat array."""
	return TableReader(
		os.path.join(sim_out, listener)).readColumn(column).flatten()


@pytest.fixture(scope='module')
def at_rest_300s():
	"""A 300 s un-stimulated run (all agonists at REST), shared by the resting
	tests: the endogenous fixed point — no IP3R firing, terminal outputs idle."""
	return _run(300, **_AT_REST)


@pytest.mark.slow
class TestRestingState:
	"""(a) The model settles to and holds its resting equilibrium."""

	def test_resting_cytosolic_ca_settles(self, at_rest_300s):
		"""Un-stimulated, cytosolic Ca²⁺ relaxes to ~100 nM and holds there.

		The fixed point is a slow relaxation, not a flat line: an early
		overshoot to ~130 nM decays back into a ~90–110 nM band by 300 s, with
		a late slope under ~0.3 nM/s. "Tolerable equilibrium" = settled near
		100 nM with bounded late drift (measured: final ~95 nM; drift over
		[200, 300] s ~+4 nM).
		"""
		ca = _col(at_rest_300s, 'CalciumTrace', 'ca_cyt_nM')
		assert 85.0 < ca[-1] < 115.0            # settled near 100 nM
		assert ca.max() < 140.0                 # bounded early overshoot (~130)
		assert abs(ca[-1] - ca[200]) < 12.0     # late drift tolerable (~4 nM)

	def test_resting_store_and_ip3_held(self, at_rest_300s):
		"""The DTS store and IP₃ hold their resting levels — no slow drain or
		run-up. Measured at 300 s: DTS ~255 µM, IP₃ ~47 nM.
		"""
		dts = _col(at_rest_300s, 'CalciumTrace', 'ca_dts_uM')
		ip3 = _col(at_rest_300s, 'CalciumTrace', 'ip3_nM')
		assert 235.0 < dts[-1] < 270.0
		assert 43.0 < ip3[-1] < 53.0


@pytest.mark.slow
class TestRestingQuiescence:
	"""(g) Resting-quiescence invariant — the property every downstream layer is
	built to preserve: with no agonist, every activation-gated output is
	*exactly* zero (the PKC×Ca coincidence gate is structurally 0 at rest). One
	guard protecting the "normalised at rest" design across secretion,
	thromboxane, integrin, the PI3K/Akt arm, and the receptor-feedback brakes.
	Exact-zero, so no band is needed — this is the most robust check in the
	model and the reason every new module leaves the Dolan benchmark unmoved."""

	def test_all_terminal_outputs_zero_at_rest(self, at_rest_300s):
		s = at_rest_300s
		assert _col(s, 'ThromboxaneTrace', 'txa2_uM').max() == 0.0
		assert _col(s, 'IntegrinTrace', 'active_frac').max() == 0.0
		assert _col(s, 'IntegrinTrace', 'akt_active').max() == 0.0
		assert _col(s, 'IntegrinTrace', 'rap1b_gtp').max() == 0.0
		assert _col(s, 'SecretionTrace', 'serotonin_released_frac').max() == 0.0
		assert _col(s, 'CalciumTrace', 'p2y12_active_frac').max() == 0.0
		assert _col(s, 'CalciumTrace', 'p2y1_desensitised_frac').max() == 0.0


@pytest.mark.slow
class TestDolanTransient:
	"""(b) The Dolan & Diamond 2014 Fig. 4 validation, expressed as behaviour
	rather than a 5/5 count: paired agonist-evoked transients with and without
	extracellular Ca²⁺. The +Ca condition peaks higher because store-operated
	entry (SOCE) refills the cytosol; under EDTA, SOCE is inactive. Config
	matches the scored Phase 3 run (length 200 s, no settle delay)."""

	def test_paired_peaks_and_soce_differential(self):
		"""+Ca and EDTA peaks sit in their Dolan-anchored bands, and the
		SOCE-dependent differential (the heart of Fig. 4) is large.

		Measured (seed 0): +Ca peak ~521 nM, EDTA peak ~296 nM, differential
		~225 nM. The peak bands are tighter than the ±30 % Dolan criterion so
		they also guard against regression, while staying within it.
		"""
		plus = _col(_run(200), 'CalciumTrace', 'ca_cyt_nM')
		edta = _col(_run(200, ca_ex_mM=0.0), 'CalciumTrace', 'ca_cyt_nM')
		assert 440.0 < plus.max() < 600.0       # +Ca peak (~521)
		assert 250.0 < edta.max() < 345.0       # EDTA peak (~296)
		# SOCE differential: the +Ca peak exceeds the EDTA peak substantially.
		assert plus.max() - edta.max() > 150.0  # measured ~225


@pytest.mark.slow
class TestMcuKnockout:
	"""(c) MCU knockout reduces the agonist-evoked cytosolic Ca²⁺ peak, matching
	the *direction* of the platelet MCU-knockout literature (Ghatge 2026; Ajanel
	2025). The detailed coupling-formula guards (and the resting-state-preserved
	check) live in ``test_validation_targets.py::TestMcuCouplingDirection``;
	this is the headline behavioural contract."""

	def test_knockout_lowers_peak_and_auc(self):
		"""Coupling on (shipped default): MCU KO lowers the evoked peak (~18 %)
		and AUC (~11 %). Coupling off (buffer-only): KO *raises* the peak — the
		sign-flip a deleted or inverted coupling would fail. Measured: WT ~530,
		KO ~435, KO-decoupled ~602 nM (60 s settle, length 300 s).
		"""
		wt = _col(_run(300, agonist_delay_s=60.0, mcu_vmax_scale=1.0),
			'CalciumTrace', 'ca_cyt_nM')
		ko = _col(_run(300, agonist_delay_s=60.0, mcu_vmax_scale=0.0),
			'CalciumTrace', 'ca_cyt_nM')
		ko_off = _col(_run(300, agonist_delay_s=60.0, mcu_vmax_scale=0.0,
			mito_coupling_gain=0.0), 'CalciumTrace', 'ca_cyt_nM')
		assert ko.max() < 0.90 * wt.max()       # KO lowers the peak (~18 %)
		assert ko.sum() < wt.sum()              # ... and the AUC (1 s steps)
		assert ko_off.max() > wt.max()          # buffer-only raises it (sign flip)


@pytest.mark.slow
class TestP2Y12Block:
	"""(d) P2Y₁₂ blockade (clopidogrel / ticagrelor / cangrelor) raises cAMP and
	lowers integrin activation. Because cytosolic Ca²⁺ is store-clamped, the
	block is *invisible* on the Ca²⁺ trace — it is read on cAMP / P2Y₁₂ occupancy
	(near-binary) and on the functional integrin (PAC-1) output, which shows a
	modest, time-dependent drop that converges as ADP clears. Default agonists
	include 10 µM ADP, which drives P2Y₁₂."""

	def test_block_raises_camp_and_lowers_integrin(self):
		"""Measured (length 200): cAMP min 0.26 → ~1.0 µM; P2Y₁₂ active 0.95 →
		exactly 0; PAC-1 active_frac at t=100 s 0.93 → 0.79.
		"""
		base = _run(200)
		blocked = _run(200, p2y12_block=1.0)
		# cAMP: P2Y₁₂ (Gi) lowers it when driven; the block restores basal cAMP.
		assert _col(base, 'CalciumTrace', 'camp_uM').min() < 0.5
		assert _col(blocked, 'CalciumTrace', 'camp_uM').min() > 0.95
		assert _col(blocked, 'CalciumTrace', 'p2y12_active_frac').max() == 0.0
		# Functional integrin: assert at the mid-rise (t=100 s), not the
		# endpoint — the effect converges as ADP clears (single-cell, clearing
		# agonist; see the test docstring).
		base_pac1 = _col(base, 'IntegrinTrace', 'active_frac')
		blk_pac1 = _col(blocked, 'IntegrinTrace', 'active_frac')
		assert base_pac1[100] > blk_pac1[100] + 0.05


@pytest.mark.slow
class TestAspirin:
	"""(e) Aspirin (COX-1 knockout, ``cox1_factor=0``) abolishes thromboxane
	synthesis. The PKC×Ca synthesis *gate* still fires (it is production that is
	lost), so the readout is the product, not the gate."""

	def test_aspirin_abolishes_thromboxane(self):
		"""Measured (length 100): TXA₂ peak 0.86 µM → exactly 0; TXB₂ (the stable
		metabolite) 38732 → 0."""
		base = _run(100)
		aspirin = _run(100, cox1_factor=0.0)
		assert 0.4 < _col(base, 'ThromboxaneTrace', 'txa2_uM').max() < 1.3
		assert _col(aspirin, 'ThromboxaneTrace', 'txa2_uM').max() == 0.0
		assert _col(aspirin, 'ThromboxaneTrace', 'txb2').max() == 0.0


if __name__ == '__main__':
	import sys
	sys.exit(pytest.main([__file__, '-v']))
