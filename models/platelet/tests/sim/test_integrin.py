"""Tests for IntegrinActivation + IntegrinTrace (v0.61 §3).

Minimal 2-state αIIbβ3 inside-out activation: a resting ⇌ active conformational
switch gated by PKC × Ca²⁺. The active fraction is the per-cell PAC-1 readout.
Covers:
  * resting quiescence — no activation when un-stimulated;
  * activation — PKC + Ca²⁺ drive a high active (PAC-1⁺) fraction;
  * conservation — the conformational switch is mass-neutral (count conserved);
  * Glanzmann / antagonist — integrin_act_scale = 0 abolishes activation.
"""

import os
import tempfile

import numpy as np
import pytest

from wholecell.io.tablereader import TableReader
from runscripts.manual.runPlateletSim import run_platelet_sim
from reconstruction.platelet.run_config import RunConfig
from models.platelet.processes.integrin_activation import akt_rap_step
from reconstruction.platelet.dataclasses.process.integrin_activation import (
	IntegrinActivation as IntegrinParams,
)

_TOTAL_INTEGRIN = 80000   # aIIbb3_resting[pl] initial count (species TSV)

# #73 — Akt → Rap1b arm rate constants, pulled from the dataclass so the unit
# tests track the real defaults (no drift if the constants change).
_P = IntegrinParams(sim_data=None)
_RAP_ARGS = dict(
	k_akt_on=_P.k_akt_on, k_akt_off=_P.k_akt_off,
	k_rap_form=_P.k_rap_form, k_rap_gap=_P.k_rap_gap, f_akt_gap=_P.f_akt_gap,
)


def _iterate_arm(gate, p2y12_frac, steps=80, f_akt_gap=None):
	"""Step ``akt_rap_step`` to steady state under constant inputs; return
	(akt, rap). ``f_akt_gap`` overrides the default coupling strength."""
	args = dict(_RAP_ARGS)
	if f_akt_gap is not None:
		args['f_akt_gap'] = f_akt_gap
	akt, rap = 0.0, 0.0
	for _ in range(steps):
		akt, rap = akt_rap_step(akt, rap, gate, p2y12_frac, 1.0, **args)
	return akt, rap


class TestAktRap1bArm:
	"""Fast unit tests for the pure PI3K/Akt → Rap1b step (#73), no sim.

	The non-vacuous mechanism guard is ``test_p2y12_sustains_rap`` — it fails if
	the Akt→GAP coupling is neutralised (verified by mutation: set f_akt_gap=0
	and it breaks). ``test_coupling_toggle_removes_the_difference`` is a
	companion *consistency* check (true by construction when f_akt_gap=0), not a
	regression guard.
	"""

	def test_rest_gives_no_akt_no_rap(self):
		"""No gate, no P2Y12 → Akt and Rap1b-GTP both decay to zero."""
		akt, rap = _iterate_arm(gate=0.0, p2y12_frac=0.0)
		assert akt == pytest.approx(0.0, abs=1e-9)
		assert rap == pytest.approx(0.0, abs=1e-9)

	def test_p2y12_sustains_rap(self):
		"""THE mechanism guard (non-vacuous): at a fixed inside-out gate, P2Y12
		occupancy (→ Akt) raises Rap1b-GTP above the no-P2Y12 (Akt-off) level —
		the sustain the lumped gate missed. Fails if f_akt_gap → 0 (coupling off)."""
		akt_on, rap_on = _iterate_arm(gate=0.4, p2y12_frac=1.0)
		akt_off, rap_off = _iterate_arm(gate=0.4, p2y12_frac=0.0)
		assert akt_on > 0.5          # P2Y12 drives Akt
		assert akt_off == pytest.approx(0.0, abs=1e-9)
		assert rap_on > rap_off      # Akt sustains Rap1b-GTP
		assert rap_on - rap_off > 0.1

	def test_coupling_toggle_removes_the_difference(self):
		"""Consistency check (not the guard — true by construction): with the
		Akt→GAP coupling off (f_akt_gap=0) the GAP term no longer depends on Akt,
		so the P2Y12-dependent Rap difference collapses to zero."""
		_, rap_on = _iterate_arm(gate=0.4, p2y12_frac=1.0, f_akt_gap=0.0)
		_, rap_off = _iterate_arm(gate=0.4, p2y12_frac=0.0, f_akt_gap=0.0)
		assert rap_on == pytest.approx(rap_off, abs=1e-9)

	def test_stays_in_unit_interval(self):
		"""Extreme inputs keep Akt and Rap1b-GTP within [0, 1]."""
		akt, rap = _iterate_arm(gate=1.0, p2y12_frac=1.0, f_akt_gap=1.0)
		assert 0.0 <= akt <= 1.0
		assert 0.0 <= rap <= 1.0


def _run(out_dir, length, **config_kwargs):
	"""Run one sim; config_kwargs are RunConfig fields (peaks, scales, …)."""
	run_config = RunConfig(ca_ex_mM=1.2, **config_kwargs)
	paths = run_platelet_sim(out_dir, length_sec=length, seed=0,
		log_to_shell=False, run_config=run_config)
	return paths['sim_out_dir']


def _ig(sim_out_dir, col):
	reader = TableReader(os.path.join(sim_out_dir, 'IntegrinTrace'))
	return reader.readColumn(col).flatten()


@pytest.mark.slow
class TestIntegrinActivation:

	def test_resting_quiescence_no_activation(self):
		"""Un-stimulated platelet: αIIbβ3 stays in the resting conformation."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 30, thrombin_peak_nM=0.0,
				adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
			assert _ig(sim_out, 'integrin_gate').max() == 0.0
			assert _ig(sim_out, 'active_frac').max() == 0.0
			assert _ig(sim_out, 'aIIbb3_active').max() == 0

	def test_activation_drives_pac1(self):
		"""Standard +Ca²⁺ transient drives αIIbβ3 to the high-affinity state."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100)
			gate = _ig(sim_out, 'integrin_gate')
			frac = _ig(sim_out, 'active_frac')
			active = _ig(sim_out, 'aIIbb3_active')
			assert gate[0] == 0.0           # zero until PKC activates
			assert frac[0] == 0.0
			assert frac.max() > 0.3         # a clear PAC-1⁺ fraction
			assert frac.max() < 1.0         # reversible — not fully saturated
			# Activation is monotonic non-decreasing while the agonist drives it.
			assert active[-1] > active[0]

	def test_total_integrin_conserved(self):
		"""resting ⇌ active is a conformational switch — count is conserved."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100)
			total = (_ig(sim_out, 'aIIbb3_active')
				+ _ig(sim_out, 'aIIbb3_resting'))
			assert np.all(total == _TOTAL_INTEGRIN)

	def test_glanzmann_knockout_abolishes_activation(self):
		"""integrin_act_scale = 0 (αIIbβ3 antagonist / Glanzmann) → no active."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100, integrin_act_scale=0.0)
			assert _ig(sim_out, 'active_frac').max() == 0.0
			assert _ig(sim_out, 'aIIbb3_active').max() == 0

	# ── #73 — PI3K/Akt → Rap1b arm (integration) ─────────────────────────

	def test_resting_no_akt_no_rap(self):
		"""Un-stimulated platelet: Akt and Rap1b-GTP stay at zero (the new
		arm preserves the resting-quiescence invariant)."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 30, thrombin_peak_nM=0.0,
				adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
			assert _ig(sim_out, 'akt_active').max() == 0.0
			assert _ig(sim_out, 'rap1b_gtp').max() == 0.0

	def test_rap1b_knockout_abolishes_activation(self):
		"""rap1b_scale = 0 (Rap1b / CalDAG-GEFI loss) → no integrin activation,
		even though the upstream gate and Akt are unaffected. Rap1b-GTP is the
		proximal driver."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100, rap1b_scale=0.0)
			assert _ig(sim_out, 'active_frac').max() == 0.0
			assert _ig(sim_out, 'aIIbb3_active').max() == 0
			# Upstream arm still runs (Rap1b forms; only the integrin step is cut).
			assert _ig(sim_out, 'rap1b_gtp').max() > 0.1

	def test_p2y12_block_slows_and_lowers_activation(self):
		"""P2Y12 blockade (cangrelor/ticagrelor) removes the Akt sustain, so
		integrin activates more slowly to a lower stimulated-phase level — the
		clinically-relevant, honest contrast (largest during the rise, not at
		the cleared-agonist endpoint).

		The active_frac / AUC contrast is *dominated by the pre-existing PKA/cAMP
		brake* (P2Y12 → Gi → ↑pka_brake, v0.7 Slice 2), not by #73 — so the
		#73-attributable assertion here is the Rap1b-GTP gap, which the PKA brake
		cannot move (it scales the integrin step, not Rap formation). The unit
		test ``test_p2y12_sustains_rap`` is the isolated mechanism guard."""
		with tempfile.TemporaryDirectory() as td:
			intact = _run(td + '/i', 100)
			blocked = _run(td + '/b', 100, p2y12_block=1.0)
			ai = _ig(intact, 'active_frac')
			ab = _ig(blocked, 'active_frac')
			# Akt is the P2Y12-driven discriminator: present intact, gone blocked.
			assert _ig(intact, 'akt_active').max() > 0.3
			assert _ig(blocked, 'akt_active').max() < 0.05
			# #73-attributable guard: the Rap1b-GTP sustain gap is purely the Akt
			# arm — the gates are ~equal, the PKA brake doesn't touch Rap — so this
			# isolates #73 from the cAMP route and collapses if f_akt_gap → 0.
			ri = _ig(intact, 'rap1b_gtp')
			rb = _ig(blocked, 'rap1b_gtp')
			assert ri[40:80].mean() > rb[40:80].mean() + 0.1
			# Both reach a clear active fraction (P2Y12 block is a brake, not a KO).
			assert ai.max() > 0.5 and ab.max() > 0.3
			# Faster rise intact: clearly higher mid-activation (t≈50).
			assert ai[50] > 1.3 * ab[50]
			# Larger stimulated-phase exposure intact (AUC over the window).
			assert np.trapz(ai) > 1.2 * np.trapz(ab)
