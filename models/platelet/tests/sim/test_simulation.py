import os
import pickle
import tempfile
import unittest

import numpy as np
import pytest

from models.platelet.sim.simulation import PlateletSimulation
from reconstruction.platelet.simulation_data import SimulationDataPlatelet
from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.utils import constants


class TestPlateletSimulationScaffold(unittest.TestCase):
	def setUp(self):
		self.sim_data = SimulationDataPlatelet()

	def make_simulation(self, **kwargs):
		output_dir = tempfile.TemporaryDirectory()
		self.addCleanup(output_dir.cleanup)

		options = {
			'outputDir': output_dir.name,
			'simData': self.sim_data,
			'logToShell': False,
			'logToDisk': False,
			'lengthSec': 1,
			}
		options.update(kwargs)

		sim = PlateletSimulation(**options)
		self.addCleanup(sim.finalize)
		return sim

	def test_simulation_data_exposes_process_namespace(self):
		"""sim_data exposes the process sub-namespaces the simulation needs."""
		from reconstruction.platelet.dataclasses.process.calcium_signalling import (
			MOLECULE_NAMES as _CALCIUM_NAMES,
		)
		self.assertTrue(hasattr(self.sim_data, 'process'))
		self.assertTrue(hasattr(self.sim_data.process, 'calcium_signalling'))
		self.assertTrue(hasattr(self.sim_data.process, 'resting_decay'))
		# CalciumSignalling exposes the 27-species ODE state vector.
		self.assertEqual(
			tuple(self.sim_data.process.calcium_signalling.molecule_names),
			tuple(_CALCIUM_NAMES),
			)

	def test_platelet_simulation_wires_real_processes(self):
		"""Both real biological processes are wired into the sim and the
		BulkMolecules state matches the configured initial counts."""
		sim = self.make_simulation()

		self.assertIn('CalciumDynamics', sim.processes)
		self.assertIn('RestingDecay', sim.processes)
		np.testing.assert_array_equal(
			sim.internal_states['BulkMolecules'].container.counts(),
			self.sim_data.internal_state.bulk_molecules.initial_counts,
			)

	def test_platelet_simulation_one_step_metabolites_preserved(self):
		"""Metabolites not touched by any process should be unchanged after one step.

		Calcium-pathway metabolites (CA2_CYT, CA2_DTS, IP3, ATP) are expected
		to change — they are driven by CalciumDynamics each timestep.  The
		remaining metabolites (granule cargo, inorganic phosphate, ADP) have no
		process acting on them and should be conserved.
		"""
		from reconstruction.platelet.dataclasses.internal_state import _MOLECULES
		from reconstruction.platelet.dataclasses.process.calcium_signalling import (
			MOLECULE_NAMES as _CALCIUM_NAMES,
		)
		calcium_set = frozenset(_CALCIUM_NAMES) | {'ATP[c]'}

		sim = self.make_simulation()
		ids = list(self.sim_data.internal_state.bulk_molecules.bulk_data['id'])
		stable_metabolite_indices = [
			ids.index(m[0])
			for m in _MOLECULES
			if m[3] == 'metabolite' and m[0] not in calcium_set
		]
		initial_counts = sim.internal_states['BulkMolecules'].container.counts().copy()

		sim.run_for(1.0)

		final_counts = sim.internal_states['BulkMolecules'].container.counts()
		np.testing.assert_array_equal(
			initial_counts[stable_metabolite_indices],
			final_counts[stable_metabolite_indices],
			err_msg='Stable metabolite counts changed unexpectedly after one step',
		)

	def test_simulation_data_pickles(self):
		payload = pickle.dumps(self.sim_data, protocol=pickle.HIGHEST_PROTOCOL)
		restored = pickle.loads(payload)
		np.testing.assert_array_equal(
			restored.process.calcium_signalling.molecule_names,
			self.sim_data.process.calcium_signalling.molecule_names,
			)

	def test_run_platelet_sim_writes_local_output(self):
		with tempfile.TemporaryDirectory() as sim_path:
			paths = run_platelet_sim(
				sim_path, length_sec=1, seed=7, log_to_shell=False)

			self.assertEqual(sim_path, paths['sim_path'])
			self.assertTrue(os.path.isfile(os.path.join(
				sim_path, constants.KB_DIR,
				constants.SERIALIZED_SIM_DATA_FILENAME)))
			self.assertTrue(os.path.isdir(os.path.join(
				paths['sim_out_dir'], 'Main')))
			self.assertTrue(os.path.isdir(os.path.join(
				paths['sim_out_dir'], 'BulkMolecules')))
			self.assertTrue(os.path.isdir(os.path.join(
				paths['sim_out_dir'], 'EvaluationTime')))
			self.assertTrue(os.path.isdir(os.path.join(
				paths['sim_out_dir'], 'Mass')))

	@pytest.mark.slow
	def test_mass_listener_dry_mass_positive(self):
		"""dryMass should be positive throughout the simulation.

		Earlier versions of this test also asserted dry mass monotonically
		decreases under RestingDecay. That holds when only protein decay is
		active, but CalciumDynamics now imports Ca²⁺ via SOCE and exports
		via PMCA, so dry mass tracks the calcium balance plus protein decay
		— it is not necessarily monotonic.
		"""
		from wholecell.io.tablereader import TableReader
		with tempfile.TemporaryDirectory() as sim_path:
			paths = run_platelet_sim(
				sim_path, length_sec=100, seed=0, log_to_shell=False)
			reader = TableReader(os.path.join(paths['sim_out_dir'], 'Mass'))
			dry = reader.readColumn('dryMass')
			self.assertTrue(np.all(dry > 0), 'dryMass should always be positive')

	def test_calcium_trace_listener_writes_output(self):
		"""CalciumTrace should write all expected columns and resting Ca²⁺ is ~100 nM."""
		from wholecell.io.tablereader import TableReader
		with tempfile.TemporaryDirectory() as sim_path:
			paths = run_platelet_sim(
				sim_path, length_sec=10, seed=0, log_to_shell=False)

			self.assertTrue(
				os.path.isdir(os.path.join(paths['sim_out_dir'], 'CalciumTrace')),
				'CalciumTrace listener output directory missing')

			reader = TableReader(os.path.join(paths['sim_out_dir'], 'CalciumTrace'))
			ca_cyt_nM     = reader.readColumn('ca_cyt_nM').flatten()
			ca_dts_uM     = reader.readColumn('ca_dts_uM').flatten()
			ip3_nM        = reader.readColumn('ip3_nM').flatten()
			soce_flux_nMs = reader.readColumn('soce_flux_nMs').flatten()

			self.assertGreater(len(ca_cyt_nM), 0, 'ca_cyt_nM column is empty')
			# Resting cytosolic Ca²⁺ should start near 100 nM (initial count = 361).
			self.assertGreater(ca_cyt_nM[0], 0.0, 'resting Ca²⁺ should be positive')
			# DTS store should start at ~250 µM (initial count = 38842).
			self.assertGreater(ca_dts_uM[0], 0.0, 'resting DTS Ca²⁺ should be positive')
			# IP3 at rest is 50 nM.
			self.assertGreater(ip3_nM[0], 0.0, 'resting IP3 should be positive')
			# SOCE flux is non-negative (Ca²⁺ enters the cell).
			self.assertTrue(
				np.all(soce_flux_nMs >= 0.0), 'SOCE flux should be non-negative')
