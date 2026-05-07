"""
CalciumDynamics process — v0.2 platelet whole-cell model.

Integrates the IP3-mediated Ca²⁺ ODE system each simulation timestep.
The ODE covers:

  * IP3R  — 6-state Markov model (Sneyd & Dufour 2002 type-2 kinetics,
             Purvis & Bhatt 2008 Table 1 rate constants)
  * SERCA — E1/E2 enzymatic cycle (Purvis 2008 / Dode 2002)
  * PMCA  — 2-state Michaelis–Menten (Caride 2007, basal kinetics)
  * SOCE  — STIM1 dimerisation + Orai1 flux (Dolan & Diamond 2014)
  * IP3   — forced time curve for v0.2 (Dolan 2014 Fig. S2); becomes a
             real upstream state variable in v0.3

The numeric ODE solver and rate constants live in:
  reconstruction/platelet/dataclasses/process/calcium_signalling.py

This process is a thin wrapper: it feeds the current molecule counts to
the solver and applies the returned integer count deltas.

IP3 forcing is active by default (ip3_forced=True). This drives the
characteristic Ca²⁺ transient using the pre-programmed IP3 time curve
from Dolan 2014 Fig. S2. In v0.3, when a real P2Y1 upstream process
produces IP3, set ip3_forced=False in the constructor.
"""

import numpy as np

import wholecell.processes.process as process
from wholecell.utils import units


class CalciumDynamics(process.Process):
	"""IP3-mediated Ca²⁺ dynamics in the resting and activated platelet."""

	_name = 'CalciumDynamics'

	# Set to False in v0.3 once the P2Y1 upstream process is in place.
	_ip3_forced = True
	# Seconds of settling time before the IP3 stimulus begins (default 0 = immediate).
	_ip3_delay = 0.0

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

	def calculateRequest(self):
		counts = self._molecules.total_counts()
		dt = self.timeStepSec()
		t_sim = self.time()

		self._delta, self._atp_cost = self._solver.molecules_to_next_time_step(
			counts, dt, t_sim, ip3_forced=self._ip3_forced,
			ip3_delay=self._ip3_delay)

		# Request molecules we need to take away (negative deltas).
		requests = np.maximum(-self._delta, 0)
		self._molecules.requestIs(requests)
		self._atp.requestIs(np.array([self._atp_cost]))

	def evolveState(self):
		self._molecules.countsInc(self._delta)
		# ATP consumed by SERCA and PMCA pumping.
		self._atp.countsDec(np.array([self._atp_cost]))
