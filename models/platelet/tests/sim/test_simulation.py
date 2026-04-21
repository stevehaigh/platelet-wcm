import os
import pickle
import tempfile
import unittest

import numpy as np

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
		self.assertTrue(hasattr(self.sim_data, 'process'))
		self.assertTrue(hasattr(self.sim_data.process, 'platelet_stub'))
		np.testing.assert_array_equal(
			self.sim_data.process.platelet_stub.molecule_names,
			self.sim_data.internal_state.bulk_molecules.bulk_data['id'],
			)

	def test_platelet_simulation_initializes_stub_process(self):
		sim = self.make_simulation()

		self.assertIn('PlateletStub', sim.processes)
		self.assertTrue(sim.processes['PlateletStub'].initialized)
		np.testing.assert_array_equal(
			sim.internal_states['BulkMolecules'].container.counts(),
			self.sim_data.internal_state.bulk_molecules.initial_counts,
			)

	def test_platelet_simulation_one_step_preserves_placeholder_counts(self):
		sim = self.make_simulation()
		initial_counts = sim.internal_states['BulkMolecules'].container.counts().copy()

		sim.run_for(1.0)

		process = sim.processes['PlateletStub']
		self.assertEqual(1, process.calculate_request_calls)
		self.assertEqual(1, process.evolve_state_calls)
		np.testing.assert_array_equal(
			initial_counts,
			sim.internal_states['BulkMolecules'].container.counts(),
			)

	def test_simulation_data_pickles(self):
		payload = pickle.dumps(self.sim_data, protocol=pickle.HIGHEST_PROTOCOL)
		restored = pickle.loads(payload)
		np.testing.assert_array_equal(
			restored.process.platelet_stub.molecule_names,
			self.sim_data.process.platelet_stub.molecule_names,
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
