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

Run conditions (extracellular Ca²⁺, agonist peaks, feedback gains,
perturbation scales) come from the per-run ``RunConfig`` on the simulation
(``sim.run_config``), read once in ``initialize`` — not from mutated module
globals.
"""

import numpy as np

import wholecell.processes.process as process
from wholecell.utils import units


class CalciumDynamics(process.Process):
	"""Receptor-driven Ca²⁺ dynamics in the resting and activated platelet."""

	_name = 'CalciumDynamics'

	def __init__(self):
		super(CalciumDynamics, self).__init__()

	def initialize(self, sim, sim_data):
		super(CalciumDynamics, self).initialize(sim, sim_data)

		self._config = sim.run_config
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

		# Per-step pericellular autocrine state, keyed by species (review #6):
		# secreted ADP[e] → P2Y1, synthesised TXA2[e] → TP. Updated discretely
		# by GranuleSecretion / ThromboxaneSynthesis.
		step_inputs = {
			'ADP[e]': float(self._adp_ex.total_counts()[0]),
			'TXA2[e]': float(self._txa2.total_counts()[0]),
		}

		self._delta, self._atp_cost = self._solver.molecules_to_next_time_step(
			counts, dt, t_sim, self._config, step_inputs)

		# Request molecules we need to take away (negative deltas).
		requests = np.maximum(-self._delta, 0)
		self._molecules.requestIs(requests)
		self._atp.requestIs(np.array([self._atp_cost]))

	def evolveState(self):
		self._molecules.countsInc(self._delta)
		# ATP consumed by SERCA and PMCA pumping.
		self._atp.countsDec(np.array([self._atp_cost]))
