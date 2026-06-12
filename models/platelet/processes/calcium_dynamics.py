"""
CalciumDynamics process — platelet whole-cell model.

Integrates the Ca²⁺ ODE system each simulation timestep. The ODE covers:

  * GPCR cascade — P2Y1 (ADP), PAR1/PAR4 (thrombin), P2X1 (ATP-gated);
                   drives a Gαq exchange/GTPase cycle → PLCβ activation
  * PI cycle    — PLCβ-catalysed PIP2 → IP3 + DAG, with lumped PIP2
                   resynthesis and IP3 / DAG degradation
  * IP3R  — Li-Rinzel 1994 reduction of de Young–Keizer 1992
             (quasi-steady m∞, one slow ODE for h)
  * SERCA — E1/E2 enzymatic cycle (Purvis 2008 / Dode 2002)
  * PMCA  — 5-state CaM-coupled (Caride 2007)
  * SOCE  — STIM1 dimerisation + Orai1 flux (Dolan & Diamond 2014)

The numeric ODE solver and rate constants live in:
  reconstruction/platelet/dataclasses/process/calcium_signalling.py

This process is a thin wrapper: it feeds the current molecule counts to
the solver and applies the returned integer count deltas.

Agonist stimulation is controlled by the three `_*_peak_*` class attrs.
Each is the peak concentration of an agonist during the activation
transient; `None` (default) uses the module-level peak constants, and
`0` yields a resting / un-stimulated sim. Override on the *class* before
constructing the simulation — mirrors the `cs_mod.CA_EX_UM` pattern.
"""

import numpy as np

import wholecell.processes.process as process
from wholecell.utils import units


class CalciumDynamics(process.Process):
	"""Receptor-driven Ca²⁺ dynamics in the resting and activated platelet."""

	_name = 'CalciumDynamics'

	# Agonist peak concentrations during the activation transient. `None`
	# uses the module-level defaults (THROMBIN_PEAK_NM, ADP_PEAK_UM,
	# ATP_EX_PEAK_UM in calcium_signalling); `0` yields a resting sim.
	_thrombin_peak_nM: float | None = None
	_adp_peak_uM:      float | None = None
	_atp_ex_peak_uM:   float | None = None
	# Seconds of settling time before the agonist stimulus begins.
	_agonist_delay = 0.0

	def __init__(self):
		super(CalciumDynamics, self).__init__()

	def initialize(self, sim, sim_data):
		super(CalciumDynamics, self).initialize(sim, sim_data)

		self._solver = sim_data.process.calcium_signalling
		self._molecules = self.bulkMoleculesView(
			self._solver.molecule_names)

		# ATP consumption is deducted each timestep.
		self._atp = self.bulkMoleculesView(
			np.array(['ATP[c]'], dtype='U30'))

		# Secreted extracellular ADP (autocrine P2Y1 drive, v0.61 Slice 2).
		# Read-only: its pericellular concentration augments the ADP forcing
		# inside the ODE. Updated discretely by GranuleSecretion.
		self._adp_ex = self.bulkMoleculesView(
			np.array(['ADP[e]'], dtype='U30'))

		# Synthesised TXA2 (autocrine TP → Gq drive, v0.61 Slice B). Read-only;
		# drives the TP receptor inside the ODE. Updated by ThromboxaneSynthesis.
		self._txa2 = self.bulkMoleculesView(
			np.array(['TXA2[e]'], dtype='U30'))

	def calculateRequest(self):
		counts = self._molecules.total_counts()
		dt = self.timeStepSec()
		t_sim = self.time()

		secreted_adp_count = float(self._adp_ex.total_counts()[0])
		secreted_txa2_count = float(self._txa2.total_counts()[0])

		self._delta, self._atp_cost = self._solver.molecules_to_next_time_step(
			counts, dt, t_sim,
			agonist_delay=self._agonist_delay,
			thrombin_peak_nM=self._thrombin_peak_nM,
			adp_peak_uM=self._adp_peak_uM,
			atp_ex_peak_uM=self._atp_ex_peak_uM,
			secreted_adp_count=secreted_adp_count,
			secreted_txa2_count=secreted_txa2_count,
		)

		# Request molecules we need to take away (negative deltas).
		requests = np.maximum(-self._delta, 0)
		self._molecules.requestIs(requests)
		self._atp.requestIs(np.array([self._atp_cost]))

	def evolveState(self):
		self._molecules.countsInc(self._delta)
		# ATP consumed by SERCA and PMCA pumping.
		self._atp.countsDec(np.array([self._atp_cost]))
