"""
Initial conditions for the platelet runtime scaffold.
"""

import reconstruction.platelet.initialization as init


def calcInitialConditions(sim, sim_data):
	"""Initialize platelet state from the minimal SimulationDataPlatelet stub."""
	assert sim._inheritedStatePath is None

	run_config = getattr(sim, 'run_config', None)
	count_overrides = run_config.count_overrides if run_config else None

	init.initialize_bulk_molecules(
		sim.internal_states['BulkMolecules'].container, sim_data,
		count_overrides=count_overrides)
	init.initialize_unique_molecules(
		sim.internal_states['UniqueMolecules'].container, sim_data)
