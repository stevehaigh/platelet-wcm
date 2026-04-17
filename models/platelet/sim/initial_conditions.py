"""
Initial conditions for the platelet runtime scaffold.
"""

import reconstruction.platelet.initialization as init


def calcInitialConditions(sim, sim_data):
	"""Initialize platelet state from the minimal SimulationDataPlatelet stub."""
	assert sim._inheritedStatePath is None

	init.initialize_bulk_molecules(
		sim.internal_states['BulkMolecules'].container, sim_data)
	init.initialize_unique_molecules(
		sim.internal_states['UniqueMolecules'].container, sim_data)
