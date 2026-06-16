"""
Run a platelet simulation from a JSON RunConfig file.

The TUI (and any caller wanting the *full* RunConfig surface, not just the
subset `runPlateletSim.py` exposes as CLI flags) writes a spec:

    {"length_sec": 60, "seed": 0, "run_config": { ...RunConfig fields... }}

and runs it with:

    PYTHONPATH=$PWD python runscripts/manual/runFromConfig.py --config run.json

A saved spec doubles as a reproducible, shareable preset. The full RunConfig is
recorded under the run's `metadata/run_config.json` for provenance.
"""

import argparse
import json
import os

import wholecell.utils.filepath as fp
from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.runPlateletSim import (
	resolve_sim_path, run_platelet_sim, write_metadata)


def load_spec(config_path):
	"""Read a run spec JSON into (length_sec, seed, RunConfig).

	The on-demand runner always streams a live CSV, so `live` is forced True.
	"""
	with open(config_path) as handle:
		spec = json.load(handle)
	length_sec = int(spec.get('length_sec', 60))
	seed = int(spec.get('seed', 0))
	rc_fields = dict(spec.get('run_config', {}))
	rc_fields['live'] = True
	return length_sec, seed, RunConfig(**rc_fields)


def build_parser():
	"""Construct the CLI parser."""
	parser = argparse.ArgumentParser(
		description='Run a platelet simulation from a JSON RunConfig.')
	parser.add_argument(
		'--config', required=True, help='Path to the run spec JSON.')
	parser.add_argument(
		'--out', default='tui_run',
		help='Output directory (relative dirs are created under out/). '
			 'Default = tui_run.')
	return parser


def main(argv=None):
	"""Run the config-file simulation CLI."""
	args = build_parser().parse_args(argv)
	length_sec, seed, run_config = load_spec(args.config)

	sim_path = resolve_sim_path(args.out)
	fp.makedirs(sim_path)
	write_metadata(
		sim_path, args.out, seed, length_sec, run_config.ca_ex_mM,
		thrombin_peak_nM=run_config.thrombin_peak_nM,
		adp_peak_uM=run_config.adp_peak_uM,
		atp_ex_peak_uM=run_config.atp_ex_peak_uM,
		agonist_delay=run_config.agonist_delay_s)

	paths = run_platelet_sim(
		sim_path, length_sec=length_sec, seed=seed,
		log_to_shell=False, run_config=run_config)

	# Record the full config (incl. gains/scales that write_metadata omits).
	meta_dir = fp.makedirs(sim_path, 'metadata')
	fp.write_json_file(
		os.path.join(meta_dir, 'run_config.json'), run_config.to_metadata())

	print('simOut: {}'.format(paths['sim_out_dir']))
	if paths['live_path']:
		print('Live CSV: {}'.format(paths['live_path']))


if __name__ == '__main__':
	main()
