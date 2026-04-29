import os
import tempfile
import unittest

from runscripts.manual.analysisPlatelet import (
	expand_plot_names, run_platelet_analysis)
from runscripts.manual.runPlateletSim import run_platelet_sim


class TestPlateletAnalysis(unittest.TestCase):
	def test_expand_plot_names_default(self):
		self.assertEqual(
			['scaffold_summary.py', 'calcium_trace.py'],
			expand_plot_names([]))

	def test_run_platelet_analysis_writes_plot_output(self):
		with tempfile.TemporaryDirectory() as sim_path:
			run_platelet_sim(sim_path, length_sec=1, seed=3, log_to_shell=False)

			result = run_platelet_analysis(sim_path)

			self.assertEqual(sim_path, result['sim_path'])
			self.assertIn('scaffold_summary', result['plots'])
			self.assertIn('calcium_trace', result['plots'])
			self.assertTrue(os.path.isdir(result['plot_out_dir']))
			for plot_name in ('scaffold_summary', 'calcium_trace'):
				self.assertTrue(os.path.isfile(os.path.join(
					result['plot_out_dir'], f'{plot_name}.pdf')),
					f'{plot_name}.pdf not found')
				self.assertTrue(os.path.isfile(os.path.join(
					result['plot_out_dir'], 'svg_plots', f'{plot_name}.svg')),
					f'{plot_name}.svg not found')
				self.assertTrue(os.path.isfile(os.path.join(
					result['plot_out_dir'], 'low_res_plots', f'{plot_name}.png')),
					f'{plot_name}.png not found')
