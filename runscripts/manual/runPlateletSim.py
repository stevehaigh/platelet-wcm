"""
Run a local platelet simulation without ParCa or Docker.

This writes a minimal `SimulationDataPlatelet` object to `kb/simData.cPickle`
and then runs `PlateletSimulation` directly into `simOut/`.
"""

import argparse
import os
import pickle
import sys

from models.platelet.listeners.calcium_trace import CalciumTrace
from models.platelet.processes.calcium_dynamics import CalciumDynamics
from models.platelet.sim.simulation import PlateletSimulation
from reconstruction.platelet.dataclasses.process import calcium_signalling as cs_mod
from reconstruction.platelet.simulation_data import SimulationDataPlatelet
from wholecell.utils import constants
import wholecell.utils.filepath as fp


DEFAULT_OUTDIR = 'platelet_manual'
DEFAULT_CA_EX_MM = 1.2          # Dolan 2014 nominal extracellular [Ca²⁺]
DEFAULT_AGONIST_DELAY = 0.0     # seconds of settling time before agonist onset (0 = immediate)
# Agonist peaks default to None on the CLI/API → cs_mod module constants
# (THROMBIN_PEAK_NM, ADP_PEAK_UM, ATP_EX_PEAK_UM) are used at call time.
# Pass --at-rest (or 0 for each) to disable a given receptor's stimulus.


def resolve_sim_path(sim_outdir):
	"""Resolve a user-provided output directory to an absolute path."""
	if os.path.isabs(sim_outdir):
		return sim_outdir
	if sim_outdir.startswith('out/'):
		return os.path.join(fp.ROOT_PATH, sim_outdir)
	return os.path.join(fp.ROOT_PATH, 'out', sim_outdir)


def write_metadata(sim_path, description, seed, length_sec, ca_ex_mM,
		thrombin_peak_nM=None, adp_peak_uM=None, atp_ex_peak_uM=None,
		agonist_delay=0.0):
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
		'thrombin_peak_nM': thrombin_peak_nM,
		'adp_peak_uM': adp_peak_uM,
		'atp_ex_peak_uM': atp_ex_peak_uM,
		'agonist_delay': agonist_delay,
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
		ca_ex_mM=DEFAULT_CA_EX_MM,
		thrombin_peak_nM=None, adp_peak_uM=None, atp_ex_peak_uM=None,
		agonist_delay=DEFAULT_AGONIST_DELAY, live=False):
	"""Create sim_data, write it to disk, and run a platelet simulation.

	Output is written to the E. coli-compatible nested structure so the
	webapp's Inspect Data tab can discover platelet runs automatically:
	  {sim_path}/platelet_stub_{seed:06d}/{seed:06d}/generation_000000/000000/simOut/

	Args:
		ca_ex_mM: extracellular Ca²⁺ concentration (Dolan 2014 default = 1.2 mM).
			Set to 0 to reproduce the Dolan Fig. 4 EDTA condition.
		thrombin_peak_nM / adp_peak_uM / atp_ex_peak_uM: peak agonist
			concentrations during the activation transient. `None`
			(default) uses the cs_mod module-level peak constants
			(THROMBIN_PEAK_NM = 1, ADP_PEAK_UM = 10, ATP_EX_PEAK_UM = 10).
			Pass `0` to disable a given receptor's stimulus (all three
			zero = a resting sim with no extracellular stimulus).
		agonist_delay: seconds before the agonist stimulus begins. The
			model settles at its fixed point during [0, agonist_delay)
			before the time courses kick in. Default 0 = immediate.
	"""
	# Override the module-level CA_EX_UM (in µM) before constructing the
	# CalciumSignalling solver. Phase 3 path A: this is how the
	# Dolan Fig. 4 with-vs-without extracellular Ca²⁺ comparison is run.
	cs_mod.CA_EX_UM = float(ca_ex_mM) * 1000.0

	# Plumb the agonist peaks and delay onto the process class before the
	# sim is constructed. Mirrors the cs_mod.CA_EX_UM override pattern.
	CalciumDynamics._thrombin_peak_nM = thrombin_peak_nM
	CalciumDynamics._adp_peak_uM = adp_peak_uM
	CalciumDynamics._atp_ex_peak_uM = atp_ex_peak_uM
	CalciumDynamics._agonist_delay = float(agonist_delay)

	# Optionally write a live CSV alongside the binary output so a viewer
	# script can tail it in real time.
	live_path = None
	if live:
		live_path = os.path.join(sim_path, 'live.csv')
		CalciumTrace._live_path = live_path
	else:
		CalciumTrace._live_path = None

	sim_data = SimulationDataPlatelet()
	write_sim_data(sim_path, sim_data)

	# Per-run metadata so downstream tooling (replayTui, analysis) can
	# recover the agonist peaks + extracellular-Ca²⁺ condition. Issue #47
	# fix: this lived in the CLI `main()` only, so dose-sweep cells (which
	# call run_platelet_sim() directly) had no metadata.json and
	# replayTui silently fell back to the module-default agonist peaks.
	description = os.path.basename(os.path.normpath(sim_path)) or sim_path
	write_metadata(sim_path, description, seed, length_sec, ca_ex_mM,
		thrombin_peak_nM=thrombin_peak_nM,
		adp_peak_uM=adp_peak_uM,
		atp_ex_peak_uM=atp_ex_peak_uM,
		agonist_delay=agonist_delay)

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
		'live_path': live_path,
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
	parser.add_argument(
		'--at-rest',
		dest='at_rest',
		action='store_true',
		help=('Shorthand: zero all three agonist peaks so the cell sits at '
			  'its endogenous fixed point with no extracellular stimulus. '
			  'Equivalent to --thrombin-peak-nM 0 --adp-peak-uM 0 '
			  '--atp-ex-peak-uM 0.'))
	parser.add_argument(
		'--thrombin-peak-nM',
		dest='thrombin_peak_nM',
		type=float,
		default=None,
		help=('Peak thrombin (nM) during the activation transient. '
			  'Default = module-level THROMBIN_PEAK_NM (1.0 nM). '
			  'Set 0 to disable PAR1/PAR4 stimulation.'))
	parser.add_argument(
		'--adp-peak-uM',
		dest='adp_peak_uM',
		type=float,
		default=None,
		help=('Peak ADP (µM) during the activation transient. '
			  'Default = module-level ADP_PEAK_UM (10 µM). '
			  'Set 0 to disable P2Y1 stimulation.'))
	parser.add_argument(
		'--atp-ex-peak-uM',
		dest='atp_ex_peak_uM',
		type=float,
		default=None,
		help=('Peak extracellular ATP (µM) during the activation transient. '
			  'Default = module-level ATP_EX_PEAK_UM (10 µM). '
			  'Set 0 to disable P2X1 stimulation.'))
	parser.add_argument(
		'--agonist-delay',
		dest='agonist_delay',
		type=float,
		default=DEFAULT_AGONIST_DELAY,
		help=('Seconds before the agonist stimulus begins. The model '
			  'settles at its fixed point during this period. Default = '
			  '{:.0f} (immediate stimulus).'.format(DEFAULT_AGONIST_DELAY)))
	parser.add_argument(
		'--live',
		dest='live',
		action='store_true',
		help=('Write a live CSV (live.csv) to the sim output directory each '
			  'timestep so a viewer script can plot the run in real time. '
			  'Use with runscripts/manual/livePlot.py.'))
	parser.set_defaults(log_to_shell=True, at_rest=False)
	return parser


def main(argv=None):
	"""Run the platelet simulation CLI."""
	args = build_parser().parse_args(argv)
	sim_outdir = args.out_override or args.sim_outdir
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)

	# --at-rest is shorthand for zeroing all three agonist peaks. Explicit
	# per-receptor overrides still win if combined (e.g. --at-rest
	# --thrombin-peak-nM 2 for a thrombin-only stimulation).
	thrombin_peak_nM = args.thrombin_peak_nM
	adp_peak_uM = args.adp_peak_uM
	atp_ex_peak_uM = args.atp_ex_peak_uM
	if args.at_rest:
		if thrombin_peak_nM is None:
			thrombin_peak_nM = 0.0
		if adp_peak_uM is None:
			adp_peak_uM = 0.0
		if atp_ex_peak_uM is None:
			atp_ex_peak_uM = 0.0

	# metadata is now written inside run_platelet_sim().
	paths = run_platelet_sim(
		sim_path,
		length_sec=args.length_sec,
		seed=args.seed,
		log_to_shell=args.log_to_shell,
		ca_ex_mM=args.ca_ex_mM,
		thrombin_peak_nM=thrombin_peak_nM,
		adp_peak_uM=adp_peak_uM,
		atp_ex_peak_uM=atp_ex_peak_uM,
		agonist_delay=args.agonist_delay,
		live=args.live,
		)

	print('Wrote platelet run to {}'.format(paths['sim_path']))
	print('simOut: {}'.format(paths['sim_out_dir']))
	if args.ca_ex_mM != DEFAULT_CA_EX_MM:
		print('Extracellular Ca²⁺ override: {:.2f} mM'.format(args.ca_ex_mM))
	all_zero = (thrombin_peak_nM == 0.0 and adp_peak_uM == 0.0
		and atp_ex_peak_uM == 0.0)
	if all_zero:
		print('All agonist peaks zero — running at rest (no stimulus)')
	else:
		for label, val, unit, default_const in (
				('Thrombin', thrombin_peak_nM, 'nM', cs_mod.THROMBIN_PEAK_NM),
				('ADP',      adp_peak_uM,      'µM', cs_mod.ADP_PEAK_UM),
				('ATP_ex',   atp_ex_peak_uM,   'µM', cs_mod.ATP_EX_PEAK_UM)):
			if val is not None and val != default_const:
				print(f'{label} peak override: {val:g} {unit} '
					f'(default {default_const:g} {unit})')
	if args.agonist_delay > 0:
		print('Agonist stimulus delayed by {:.0f} s'.format(args.agonist_delay))
	if args.live and paths['live_path']:
		print('Live CSV: {}'.format(paths['live_path']))
		print('Viewer:   PYTHONPATH=$PWD pyenv exec python '
			  'runscripts/manual/livePlot.py {}'.format(paths['live_path']))


if __name__ == '__main__':
	main()
