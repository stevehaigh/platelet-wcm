"""
Mass listener for the platelet whole-cell model.

Records total dry mass, protein mass, and metabolite mass each timestep.
Mass units are femtograms (fg) throughout — the same unit used in BulkMolecules.

Under resting conditions the only active process is RestingDecay, so
`dryMass` should decrease monotonically as proteins turn over.
"""

import numpy as np

import wholecell.listeners.listener


class Mass(wholecell.listeners.listener.Listener):
	"""Record cellular mass by category each timestep."""

	_name = 'Mass'

	def __init__(self, *args, **kwargs):
		self.internal_states = None
		self.massUnits = 'fg'
		super().__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super().initialize(sim, sim_data)

		self.internal_states = sim.internal_states

		self.proteinIndex = sim_data.submass_name_to_index['protein']
		self.metaboliteIndex = sim_data.submass_name_to_index['metabolite']

		self.dryMass = 0.0
		self.proteinMass = 0.0
		self.metaboliteMass = 0.0
		self.growth = np.nan
		self._setInitial = False

		self.registerLoggedQuantity('Dry mass\n(fg)', 'dryMass', '.1f')
		self.registerLoggedQuantity('Protein\n(fg)', 'proteinMass', '.1f')
		self.registerLoggedQuantity('Growth\n(fg/s)', 'growth', '.4f')

	def update(self):
		old_dry = self.dryMass

		all_submasses = sum(
			state.mass() for state in self.internal_states.values())

		self.proteinMass = float(all_submasses[self.proteinIndex])
		self.metaboliteMass = float(all_submasses[self.metaboliteIndex])
		self.dryMass = self.proteinMass + self.metaboliteMass

		if not self._setInitial:
			self._setInitial = True
			self.dryMassInitial = self.dryMass
			self.proteinMassInitial = self.proteinMass
			self.metaboliteMassInitial = self.metaboliteMass

		self.growth = self.dryMass - old_dry if self.simulationStep() > 0 else np.nan
		self.dryMassFoldChange = self.dryMass / self.dryMassInitial

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			mass_units=self.massUnits,
			submass_names=['protein', 'metabolite'],
		)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time=self.time(),
			simulationStep=self.simulationStep(),
			dryMass=self.dryMass,
			proteinMass=self.proteinMass,
			metaboliteMass=self.metaboliteMass,
			growth=self.growth,
			dryMassFoldChange=self.dryMassFoldChange,
		)
