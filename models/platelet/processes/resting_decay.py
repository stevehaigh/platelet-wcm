"""
Resting-state protein decay process for the platelet whole-cell model.

Anucleate platelets cannot replenish proteins; all proteins decay exponentially
with a half-life of ~7 days (Burkhart 2012).

Each timestep, each molecule decays independently with probability:
    p = 1 - exp(-ln2 * dt / t_half)

Losses are drawn from Binomial(count, p) — correct for integer molecule counts.
"""

import numpy as np

import wholecell.processes.process


class RestingDecay(wholecell.processes.process.Process):
	"""Exponential protein decay in the resting (anucleate) platelet."""

	_name = 'RestingDecay'

	def __init__(self):
		super(RestingDecay, self).__init__()

	def initialize(self, sim, sim_data):
		super(RestingDecay, self).initialize(sim, sim_data)
		dc = sim_data.process.resting_decay
		self._half_life = dc.protein_half_life
		self._molecules = self.bulkMoleculesView(dc.molecule_names)

	def calculateRequest(self):
		self._molecules.requestAll()

	def evolveState(self):
		counts = self._molecules.counts()
		dt = self.timeStepSec()
		decay_prob = 1.0 - np.exp(-np.log(2) * dt / self._half_life)
		lost = np.random.binomial(counts, decay_prob)
		self._molecules.countsDec(lost)
