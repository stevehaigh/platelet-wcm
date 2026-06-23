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
from dataclasses import dataclass, field
from typing import Dict


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

	p2y12_block: float = 0.0
	"""P2Y12 competitive-antagonist occupancy (v0.7 Slice 2, issue #10);
	0 = no drug (P2Y12 fully active), 1 = full block (cangrelor / ticagrelor /
	clopidogrel). Scales the ADP on-rate at P2Y12, so blockade keeps cAMP at
	its basal tone → PKA brake intact → reduced activation."""

	# ── cAMP-raising (Gs) arm — antiplatelet/vasodilator drugs (v0.7 Slice 1) ─
	pgi2_nM: float = 0.0
	"""PGI₂ / iloprost concentration (nM) — Gs → adenylyl cyclase stimulation
	that RAISES cAMP (v0.7 Slice 1). 0 = none. The prostacyclin arm; suppresses
	activation (↑PKA brake → ↓integrin/secretion, ↑VASP-P)."""

	forskolin: float = 0.0
	"""Forskolin — direct adenylyl-cyclase activator, as a fold-increase in AC
	production (v0.7 Slice 1). 0 = none. Tool compound; raises cAMP."""

	pde3_block: float = 0.0
	"""PDE3A inhibition fraction (cilostazol / dipyridamole) ∈ [0,1] (v0.7
	Slice 1). 0 = none; reduces cAMP degradation → raises cAMP."""

	# ── Single-constant perturbation scales (1.0 = baseline; 0.0 = knockout) ─
	k_des_scale: float = 1.0
	"""Scale on PKC → P2Y1 desensitisation rate (K_P2Y1_DES['k_des'])."""

	k_plcb_phos_scale: float = 1.0
	"""Scale on PKC → PLCβ phosphorylation rate (K_PLCB_PHOS['k_plcb_phos'])."""

	pmca_kcat_scale: float = 1.0
	"""Scale on PMCA turnover (K_PMCA['k_cat']) — the runPerturbation pmca sweep."""

	mcu_vmax_scale: float = 1.0
	"""Scale on MCU V_max (K_MITO['V_max_MCU']) — the runPerturbation mcu sweep."""

	mito_coupling_gain: float = 1.0
	"""#76 Part 2 — toggle/strength of the MCU → IP3R-relief coupling, via
	`ip3r_relief_factor`. Scales the evoked relief of the IP3R's Ca²⁺-dependent
	inactivation that MCU provides at mitochondria–DTS (MAM) contacts. 1.0 = full
	coupling; 0.0 = decoupled (Part-1-only). MCU knockout (mcu_vmax_scale→0) → relief
	lost → reduced agonist-evoked release → cytosolic Ca²⁺ reduced; the fuller store
	then lowers SOCE indirectly (SOCE is NOT gated directly). Ghatge/Ajanel 2025."""

	# ── Initial-count overrides (Tier 2 — expression knockouts) ───────────
	count_overrides: Dict[str, int] = field(default_factory=dict)
	"""Per-run overrides of resting copy numbers, keyed by molecule id (e.g.
	``{'P2Y1_inactive[pl]': 0}``). Applied when the bulk-molecule state is
	seeded (``initialize_bulk_molecules``), over the sim_data baseline. Empty
	(default) → byte-identical to the unmodified run. Setting an entity's
	copy number to 0 is an expression knockout; the TUI expands a logical
	entity (e.g. "P2Y1") into all its conformational sub-states via the
	knockout entity map before populating this."""

	def to_metadata(self) -> dict:
		"""Flat dict of the config for the run ``metadata/`` JSON (provenance)."""
		return dataclasses.asdict(self)
