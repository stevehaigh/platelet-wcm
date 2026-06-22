import os
import tempfile
import unittest

import numpy as np
import pytest

from models.platelet.analysis.single.granule_secretion import _onset_time
from runscripts.manual.analysisPlatelet import (
	expand_plot_names, run_platelet_analysis)
from runscripts.manual.runPlateletSim import run_platelet_sim


class TestPlateletAnalysis(unittest.TestCase):
	def test_expand_plot_names_default(self):
		self.assertEqual(
			['scaffold_summary.py', 'calcium_trace.py', 'granule_secretion.py'],
			expand_plot_names([]))

	@pytest.mark.slow
	def test_run_platelet_analysis_writes_plot_output(self):
		with tempfile.TemporaryDirectory() as sim_path:
			run_platelet_sim(sim_path, length_sec=1, seed=3, log_to_shell=False)

			result = run_platelet_analysis(sim_path)

			self.assertEqual(sim_path, result['sim_path'])
			self.assertIn('scaffold_summary', result['plots'])
			self.assertIn('calcium_trace', result['plots'])
			self.assertTrue(os.path.isdir(result['plot_out_dir']))
			for plot_name in ('scaffold_summary', 'calcium_trace',
					'granule_secretion'):
				self.assertTrue(os.path.isfile(os.path.join(
					result['plot_out_dir'], f'{plot_name}.pdf')),
					f'{plot_name}.pdf not found')
				self.assertTrue(os.path.isfile(os.path.join(
					result['plot_out_dir'], 'svg_plots', f'{plot_name}.svg')),
					f'{plot_name}.svg not found')
				self.assertTrue(os.path.isfile(os.path.join(
					result['plot_out_dir'], 'low_res_plots', f'{plot_name}.png')),
					f'{plot_name}.png not found')

	@pytest.mark.slow
	def test_out_name_overrides_single_plot_filename(self):
		with tempfile.TemporaryDirectory() as sim_path:
			run_platelet_sim(sim_path, length_sec=1, seed=5, log_to_shell=False)

			result = run_platelet_analysis(
				sim_path, ['demo_calcium'], out_name='custom_fig')

			# figure is written under the custom base name, not the module name
			self.assertEqual(['custom_fig'], result['plots'])
			self.assertTrue(os.path.isfile(os.path.join(
				result['plot_out_dir'], 'custom_fig.pdf')))
			self.assertFalse(os.path.isfile(os.path.join(
				result['plot_out_dir'], 'demo_calcium.pdf')))

	def test_granule_secretion_plot_in_active(self):
		self.assertIn('granule_secretion.py', expand_plot_names(['ACTIVE']))

	def test_onset_time_threshold_crossing(self):
		# Rises through the 1% onset threshold at t=2 s.
		t = np.array([0.0, 1.0, 2.0, 3.0])
		frac = np.array([0.0, 0.005, 0.05, 0.5])
		self.assertEqual(2.0, _onset_time(t, frac))

	def test_onset_time_no_release_returns_none(self):
		t = np.array([0.0, 1.0, 2.0])
		frac = np.array([0.0, 0.0, 0.0])
		self.assertIsNone(_onset_time(t, frac))
