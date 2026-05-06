"""
Run a local platelet simulation without ParCa or Docker.

This writes a minimal `SimulationDataPlatelet` object to `kb/simData.cPickle`
and then runs `PlateletSimulation` directly into `simOut/`.
"""

import argparse
import os
import pickle
import sys

from models.platelet.sim.simulation import PlateletSimulation
from reconstruction.platelet.dataclasses.process import calcium_signalling as cs_mod
from reconstruction.platelet.simulation_data import SimulationDataPlatelet
from wholecell.utils import constants
import wholecell.utils.filepath as fp


DEFAULT_OUTDIR = 'platelet_manual'
DEFAULT_CA_EX_MM = 1.2     # Dolan 2014 nominal extracellular [Ca²⁺]


def resolve_sim_path(sim_outdir):
	"""Resolve a user-provided output directory to an absolute path."""
	if os.path.isabs(sim_outdir):
		return sim_outdir
	if sim_outdir.startswith('out/'):
		return os.path.join(fp.ROOT_PATH, sim_outdir)
	return os.path.join(fp.ROOT_PATH, 'out', sim_outdir)


def write_metadata(sim_path, description, seed, length_sec, ca_ex_mM):
	"""Write a small metadata file for a local platelet run."""
	metadata = {
		'git_hash': fp.git_hash(),
		'git_branch': fp.git_branch(),
		'description': description,
		'time': fp.timestamp(),
		'python': sys.version.splitlines()[0],
		'seed': seed,
		'length_sec': length_sec,
		'ca_ex_mM': ca_ex_mM,
		'variant': 'platelet_stub',
		'analysis_type': None,
		}
	metadata_dir = fp.makedirs(sim_path, constants.METADATA_DIR)
	metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)
	fp.write_json_file(metadata_path, metadata)


def write_sim_data(sim_path, sim_data):
	"""Serialize the platelet sim_data object into the run directory."""
	kb_directory = fp.makedirs(sim_path, constants.KB_DIR)
	sim_data_path = os.path.join(
		kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
	with open(sim_data_path, 'wb') as sim_data_file:
		pickle.dump(sim_data, sim_data_file, protocol=pickle.HIGHEST_PROTOCOL)
	return sim_data_path


def run_platelet_sim(sim_path, length_sec, seed, log_to_shell=True,
		ca_ex_mM=DEFAULT_CA_EX_MM):
	"""Create sim_data, write it to disk, and run a platelet simulation.

	Output is written to the E. coli-compatible nested structure so the
	webapp's Inspect Data tab can discover platelet runs automatically:
	  {sim_path}/platelet_stub_{seed:06d}/{seed:06d}/generation_000000/000000/simOut/

	`ca_ex_mM` overrides the extracellular Ca²⁺ concentration used by the
	calcium signalling ODE (Dolan 2014 default = 1.2 mM). Setting this to
	0 reproduces the Dolan Fig. 4 EDTA condition.
	"""
	# Override the module-level CA_EX_UM (in µM) before constructing the
	# CalciumSignalling solver. Phase 3 path A: this is how the
	# Dolan Fig. 4 with-vs-without extracellular Ca²⁺ comparison is run.
	cs_mod.CA_EX_UM = float(ca_ex_mM) * 1000.0

	sim_data = SimulationDataPlatelet()
	write_sim_data(sim_path, sim_data)

	cell_dir = fp.makedirs(
		sim_path,
		f'platelet_stub_{seed:06d}',
		f'{seed:06d}',
		'generation_000000',
		'000000',
	)
	sim_out_dir = fp.makedirs(cell_dir, 'simOut')
	sim = PlateletSimulation(
		simData=sim_data,
		outputDir=sim_out_dir,
		lengthSec=length_sec,
		seed=seed,
		logToShell=log_to_shell,
		logToDisk=True,
		)
	sim.run()

	return {
		'sim_path': sim_path,
		'sim_out_dir': sim_out_dir,
		}


def build_parser():
	"""Construct the CLI parser for the platelet local runner."""
	parser = argparse.ArgumentParser(
		description='Run a local platelet scaffold simulation.')
	parser.add_argument(
		'sim_outdir',
		nargs='?',
		default=DEFAULT_OUTDIR,
		help=('Output directory for the run. Relative paths are created under '
			  '"out/". Default = {!r}.'.format(DEFAULT_OUTDIR)))
	parser.add_argument(
		'--out',
		dest='out_override',
		help='Optional explicit output directory override.')
	parser.add_argument(
		'--length', '--length-sec',
		dest='length_sec',
		type=int,
		default=60,
		help='Simulation length in seconds. Default = 60.')
	parser.add_argument(
		'--seed',
		type=int,
		default=0,
		help='Random seed for the simulation. Default = 0.')
	parser.add_argument(
		'--no-log-to-shell',
		dest='log_to_shell',
		action='store_false',
		help='Disable shell logging during the run.')
	parser.add_argument(
		'--ca-ex-mM',
		dest='ca_ex_mM',
		type=float,
		default=DEFAULT_CA_EX_MM,
		help=('Extracellular Ca²⁺ concentration in mM. '
			  'Default = {:.1f} (Dolan 2014 nominal). '
			  'Set to 0 for the Dolan Fig. 4 EDTA / no-extracellular-Ca '
			  'condition.'.format(DEFAULT_CA_EX_MM)))
	parser.set_defaults(log_to_shell=True)
	return parser


def main(argv=None):
	"""Run the platelet simulation CLI."""
	args = build_parser().parse_args(argv)
	sim_outdir = args.out_override or args.sim_outdir
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)

	description = os.path.basename(os.path.normpath(sim_path)) or sim_path
	write_metadata(sim_path, description, args.seed, args.length_sec,
		args.ca_ex_mM)
	paths = run_platelet_sim(
		sim_path,
		length_sec=args.length_sec,
		seed=args.seed,
		log_to_shell=args.log_to_shell,
		ca_ex_mM=args.ca_ex_mM,
		)

	print('Wrote platelet run to {}'.format(paths['sim_path']))
	print('simOut: {}'.format(paths['sim_out_dir']))
	if args.ca_ex_mM != DEFAULT_CA_EX_MM:
		print('Extracellular Ca²⁺ override: {:.2f} mM'.format(args.ca_ex_mM))


if __name__ == '__main__':
	main()
