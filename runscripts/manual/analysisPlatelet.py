"""Run platelet-specific analysis plots for a local platelet scaffold run."""

import importlib
import os

from models.platelet.analysis.single import TAGS
from wholecell.utils import constants, scriptBase
import wholecell.utils.filepath as fp
from wholecell.io import run_results as _results


DEFAULT_PLOTS = ['DEFAULT']


def resolve_sim_path(sim_dir):
	"""Resolve a platelet run directory to an absolute path."""
	if os.path.isabs(sim_dir):
		return sim_dir
	if sim_dir.startswith('out/'):
		return os.path.join(fp.ROOT_PATH, sim_dir)
	return os.path.join(fp.ROOT_PATH, 'out', sim_dir)


def expand_plot_names(plot_names):
	"""Expand tag names like DEFAULT or CORE into concrete plot module names."""
	expanded = []
	for plot_name in (plot_names or DEFAULT_PLOTS):
		key = plot_name.upper()
		if key in TAGS:
			expanded.extend(TAGS[key])
		else:
			name = plot_name if plot_name.endswith('.py') else plot_name + '.py'
			expanded.append(name)
	return expanded


def _find_simout(sim_path):
	"""Locate the simOut directory for a platelet run.

	Checks for the E. coli-compatible nested structure first (written by
	runPlateletSim.py), then falls back to the legacy flat structure.
	"""
	for variant in _results.find_variants(sim_path):
		cells = _results.find_cells(sim_path, variant)
		if cells:
			return cells[0]['simout_path']
	return os.path.join(sim_path, 'simOut')


def run_platelet_analysis(sim_path, plot_names=None, out_name=None):
	"""Run one or more platelet analysis plots for a local run.

	``out_name`` overrides the output file's base name; it only applies when a
	single plot is requested (otherwise plots would overwrite each other), so
	the figure can be named by e.g. a TUI save-as name instead of the module.
	"""
	sim_out_dir = _find_simout(sim_path)
	# Write plots alongside simOut/ so the web app can find them via find_plot_images
	plot_out_dir = fp.makedirs(sim_out_dir.replace('simOut', constants.PLOTOUT_DIR))
	sim_data_file = os.path.join(
		sim_path, constants.KB_DIR, constants.SERIALIZED_SIM_DATA_FILENAME)
	metadata_path = os.path.join(
		sim_path, constants.METADATA_DIR, constants.JSON_METADATA_FILE)
	metadata = fp.read_json_file(metadata_path) if os.path.exists(metadata_path) else {}
	if metadata.get('analysis_type') is None:
		metadata.pop('analysis_type', None)

	expanded = expand_plot_names(plot_names)
	plots = []
	for filename in expanded:
		module_name = os.path.splitext(filename)[0]
		# Custom name only when a single plot is requested, else plots collide.
		file_name = out_name if (out_name and len(expanded) == 1) else module_name
		module = importlib.import_module(
			'models.platelet.analysis.single.' + module_name)
		module.Plot.main(
			sim_out_dir,
			plot_out_dir,
			file_name,
			sim_data_file,
			None,
			metadata,
			)
		plots.append(file_name)

	return {
		'sim_path': sim_path,
		'plot_out_dir': plot_out_dir,
		'plots': plots,
		}


class AnalysisPlatelet(scriptBase.ScriptBase):
	"""Run platelet analysis plots for a local platelet scaffold run."""

	def description(self):
		return 'Platelet scaffold analysis'

	def define_parameters(self, parser):
		super(AnalysisPlatelet, self).define_parameters(parser)
		parser.add_argument(
			'sim_dir',
			nargs='?',
			default='platelet_manual',
			help='Platelet run directory. Relative paths are resolved under out/.')
		parser.add_argument(
			'-p', '--plot',
			nargs='+',
			default=[],
			help='Platelet plot names or tags (DEFAULT, CORE, ACTIVE).')
		parser.add_argument(
			'--out-name',
			default=None,
			help='Override the output figure base name (single --plot only).')

	def parse_args(self):
		args = super(AnalysisPlatelet, self).parse_args()
		args.sim_path = resolve_sim_path(args.sim_dir)
		return args

	def run(self, args):
		result = run_platelet_analysis(
			args.sim_path, args.plot, out_name=getattr(args, 'out_name', None))
		print('Wrote platelet plots to {}'.format(result['plot_out_dir']))
		print('Plots: {}'.format(', '.join(result['plots'])))


if __name__ == '__main__':
	AnalysisPlatelet().cli()
