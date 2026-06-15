"""
RunConfig — per-simulation run conditions for the platelet whole-cell model.

This is **run configuration**, deliberately kept separate from the
reconstruction (``SimulationDataPlatelet`` / ``sim_data``). ``sim_data`` is the
fixed biology — rate constants, species, volumes — built once and reused across
many runs (the analogue of the ParCa output upstream). ``RunConfig`` is "which
experiment am I running": extracellular Ca²⁺, agonist doses, feedback-loop
gains, and single-constant perturbation scales.

It replaces the previous mechanism of reassigning module-level globals in place
(``cs_mod.CA_EX_UM = …``, ``CalciumDynamics._adp_peak_uM = …``,
``tx_mod.COX1_FACTOR = …``, ``K_P2Y1_DES['k_des'] = …``), which was
process-global (two simulations in one process would collide), correct only by
save/restore discipline, and invisible at the call site. A ``RunConfig`` is
constructed once, attached to the simulation, read by the processes/ODE, and
recorded in the run ``metadata/`` for provenance.

Frozen (immutable): a run's configuration cannot change mid-run, and two runs in
the same process get two independent configs.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass


@dataclass(frozen=True)
class RunConfig:
	"""Immutable per-run configuration for a platelet simulation.

	Agonist peaks use ``None`` to mean "use the module default loaded from the
	kinetics TOML" (resolved where the agonist time courses are evaluated); pass
	``0.0`` for a receptor that sees only its REST level. The ``*_scale`` fields
	multiply a single rate constant (``1.0`` = baseline, ``0.0`` = knockout) and
	replace ``runPerturbation``'s in-place mutation of the ``K_*`` dicts.
	"""

	# ── Run conditions ────────────────────────────────────────────────────
	ca_ex_mM: float = 1.2
	"""Extracellular Ca²⁺ (mM). 0 → the Dolan Fig. 4 EDTA condition (SOCE and
	the plasma-membrane leak are gated off)."""

	thrombin_peak_nM: float | None = None
	"""Peak thrombin (nM) during the activation transient; None → TOML default,
	0 → REST (PAR1/PAR4 unstimulated)."""

	adp_peak_uM: float | None = None
	"""Peak ADP (µM); None → TOML default, 0 → REST (P2Y1 unstimulated)."""

	atp_ex_peak_uM: float | None = None
	"""Peak extracellular ATP (µM); None → TOML default, 0 → REST (P2X1)."""

	agonist_delay_s: float = 0.0
	"""Seconds the model settles at its fixed point before the agonist time
	courses begin."""

	live: bool = False
	"""Write a live CSV alongside the binary output (for the live-plot viewer)."""

	# ── Feedback-loop gains (1.0 = on; 0.0 = open loop) ───────────────────
	autocrine_adp_gain: float = 1.0
	"""Gain on the autocrine secreted-ADP → P2Y1 loop (v0.61 Slice 2)."""

	autocrine_txa2_gain: float = 1.0
	"""Gain on the autocrine TXA₂ → TP → Gq loop (v0.61 Slice B)."""

	cox1_factor: float = 1.0
	"""COX-1 availability for TXA₂ synthesis — the aspirin knob (0 = aspirin)."""

	integrin_act_scale: float = 1.0
	"""Scale on αIIbβ3 inside-out activation rate (v0.61 §3); 0 = αIIbβ3
	antagonist / Glanzmann thrombasthenia (no high-affinity integrin)."""

	# ── Single-constant perturbation scales (1.0 = baseline; 0.0 = knockout) ─
	k_des_scale: float = 1.0
	"""Scale on PKC → P2Y1 desensitisation rate (K_P2Y1_DES['k_des'])."""

	k_plcb_phos_scale: float = 1.0
	"""Scale on PKC → PLCβ phosphorylation rate (K_PLCB_PHOS['k_plcb_phos'])."""

	pmca_kcat_scale: float = 1.0
	"""Scale on PMCA turnover (K_PMCA['k_cat']) — the runPerturbation pmca sweep."""

	mcu_vmax_scale: float = 1.0
	"""Scale on MCU V_max (K_MITO['V_max_MCU']) — the runPerturbation mcu sweep."""

	def to_metadata(self) -> dict:
		"""Flat dict of the config for the run ``metadata/`` JSON (provenance)."""
		return dataclasses.asdict(self)
