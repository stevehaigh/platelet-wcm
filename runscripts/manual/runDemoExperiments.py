"""Regenerate the four demo-experiment figure sets in one pass.

The demo experiments (`reports/experiments/{1..4}-*.qmd`) each contrast a
perturbation against a baseline run and embed a figure under
`reports/figures/demos/`. This driver is the single source of truth for those
runs: it declares each run's `RunConfig`, runs the simulations once, renders the
demo figures (wiring the grey baseline overlay via `PLATELET_BASELINE_SIMOUT`),
copies the PNGs into `reports/figures/demos/`, and prints a metrics table for
the prose in the `.qmd` write-ups.

Every run uses the same protocol: a `SETTLE_S`-second settle at the resting
fixed point with all agonists at REST, then the agonists switch on and the sim
runs for a further `RUN_S` seconds (so `length = SETTLE_S + RUN_S`).

    PYTHONPATH=$PWD python runscripts/manual/runDemoExperiments.py
    PYTHONPATH=$PWD python runscripts/manual/runDemoExperiments.py --only mcu_calcium

Runs land under `out/demo_experiments/<run>/`; figures are copied to
`reports/figures/demos/<figure>.png`.
"""

import argparse
import os
import shutil

import wholecell.utils.filepath as fp
from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.analysisPlatelet import run_platelet_analysis
from runscripts.manual.runPlateletSim import (
	resolve_sim_path, run_platelet_sim, write_metadata)
from wholecell.analysis.analysis_tools import LOW_RES_DIR
from wholecell.io.tablereader import TableReader

# --- experiment protocol: settle at rest, then stimulate -------------------
SETTLE_S = 300          # seconds held at the resting fixed point (no agonist)
RUN_S = 300             # seconds of agonist stimulation after the settle
LENGTH_S = SETTLE_S + RUN_S
SEED = 0

DEMOS_DIR = os.path.join(fp.ROOT_PATH, 'reports', 'figures', 'demos')
OUT_PREFIX = 'demo_experiments'

# --- the distinct simulations (run once each, keyed by name) ---------------
# `agonist_delay_s` is set to SETTLE_S for every stimulated run; the resting run
# leaves all agonists at REST so the delay is irrelevant. Omitted agonist peaks
# fall back to the Dolan standard (thrombin 1 nM, ADP 10 µM, ATP 10 µM).
RUNS = {
	'rest':  dict(thrombin_peak_nM=0, adp_peak_uM=0, atp_ex_peak_uM=0),
	'base':  dict(agonist_delay_s=SETTLE_S),
	'wbase': dict(agonist_delay_s=SETTLE_S,
				  thrombin_peak_nM=0, adp_peak_uM=0.5, atp_ex_peak_uM=0),
	'asp':   dict(agonist_delay_s=SETTLE_S, cox1_factor=0.0),
	'mcu':   dict(agonist_delay_s=SETTLE_S, mcu_vmax_scale=0.0),
	'clop':  dict(agonist_delay_s=SETTLE_S, p2y12_block=1.0),
	'wclop': dict(agonist_delay_s=SETTLE_S,
				  thrombin_peak_nM=0, adp_peak_uM=0.5, atp_ex_peak_uM=0,
				  p2y12_block=1.0),
}

# --- figure -> (run, plot module, optional baseline run for the grey overlay)
FIGURES = [
	('resting_calcium',             'rest',  'demo_calcium',     None),
	('baseline_calcium',            'base',  'demo_calcium',     None),
	('baseline_integrin',           'base',  'demo_integrin',    None),
	('aspirin_thromboxane',         'asp',   'demo_thromboxane', 'base'),
	('mcu_calcium',                 'mcu',   'demo_calcium',     'base'),
	('clopidogrel_strong_integrin', 'clop',  'demo_integrin',    'base'),
	('clopidogrel_weak_integrin',   'wclop', 'demo_integrin',    'wbase'),
]


def run_one(run_key):
	"""Run a single named simulation; return its simOut directory."""
	run_config = RunConfig(live=False, **RUNS[run_key])
	sim_path = resolve_sim_path(os.path.join(OUT_PREFIX, run_key))
	fp.makedirs(sim_path)
	write_metadata(
		sim_path, '{}/{}'.format(OUT_PREFIX, run_key), SEED, LENGTH_S,
		run_config.ca_ex_mM,
		thrombin_peak_nM=run_config.thrombin_peak_nM,
		adp_peak_uM=run_config.adp_peak_uM,
		atp_ex_peak_uM=run_config.atp_ex_peak_uM,
		agonist_delay=run_config.agonist_delay_s)
	paths = run_platelet_sim(
		sim_path, length_sec=LENGTH_S, seed=SEED,
		log_to_shell=False, run_config=run_config)
	meta_dir = fp.makedirs(sim_path, 'metadata')
	fp.write_json_file(
		os.path.join(meta_dir, 'run_config.json'), run_config.to_metadata())
	return paths['sim_out_dir']


def render_figure(name, run_key, plot_module, baseline_run, simouts):
	"""Render one demo figure and copy its PNG into reports/figures/demos/."""
	prev = os.environ.pop('PLATELET_BASELINE_SIMOUT', None)
	if baseline_run is not None:
		os.environ['PLATELET_BASELINE_SIMOUT'] = simouts[baseline_run]
	try:
		result = run_platelet_analysis(
			resolve_sim_path(os.path.join(OUT_PREFIX, run_key)),
			[plot_module], out_name=name)
	finally:
		os.environ.pop('PLATELET_BASELINE_SIMOUT', None)
		if prev is not None:
			os.environ['PLATELET_BASELINE_SIMOUT'] = prev
	src = os.path.join(result['plot_out_dir'], LOW_RES_DIR, name + '.png')
	dst = os.path.join(DEMOS_DIR, name + '.png')
	shutil.copyfile(src, dst)
	return dst


def _col(simout, listener, column):
	"""Read one listener column, or None if absent."""
	try:
		return TableReader(os.path.join(simout, listener)).readColumn(column)
	except (IOError, OSError, ValueError, KeyError):
		return None


def print_metrics(simouts):
	"""Print a compact per-run metrics table for the write-up prose."""
	hdr = ('run', 'Ca peak nM', 'Ca final nM', 'IP3 final nM',
		   'PAC-1 final', 'TXA2 peak uM')
	print('\n{:<7} {:>11} {:>12} {:>13} {:>12} {:>13}'.format(*hdr))
	for key in RUNS:
		so = simouts.get(key)
		if so is None:
			continue
		ca = _col(so, 'CalciumTrace', 'ca_cyt_nM')
		ip3 = _col(so, 'CalciumTrace', 'ip3_nM')
		paf = _col(so, 'IntegrinTrace', 'active_frac')
		txa2 = _col(so, 'ThromboxaneTrace', 'txa2_uM')
		def f(v, fmt):
			return fmt.format(v) if v is not None else '—'
		print('{:<7} {:>11} {:>12} {:>13} {:>12} {:>13}'.format(
			key,
			f(ca.max() if ca is not None else None, '{:.1f}'),
			f(ca[-1] if ca is not None else None, '{:.1f}'),
			f(ip3[-1] if ip3 is not None else None, '{:.1f}'),
			f(paf[-1] if paf is not None else None, '{:.3f}'),
			f(txa2.max() if txa2 is not None else None, '{:.3f}')))


def build_parser():
	"""Construct the CLI parser."""
	parser = argparse.ArgumentParser(
		description='Regenerate the demo-experiment figures (settle then '
					'stimulate; settle={}s, run={}s).'.format(SETTLE_S, RUN_S))
	parser.add_argument(
		'--only', nargs='+', metavar='FIGURE',
		choices=[f[0] for f in FIGURES],
		help='Render only these figures (still runs the sims they need).')
	parser.add_argument(
		'--no-metrics', action='store_true',
		help='Skip the metrics summary table.')
	return parser


def main(argv=None):
	"""Run the requested demo experiments and regenerate their figures."""
	args = build_parser().parse_args(argv)
	figures = [f for f in FIGURES if not args.only or f[0] in args.only]
	needed = {run for _, run, _, _ in figures}
	needed |= {bl for _, _, _, bl in figures if bl is not None}

	simouts = {}
	for key in RUNS:                      # deterministic order
		if key in needed:
			print('>>> sim {}'.format(key))
			simouts[key] = run_one(key)

	for name, run_key, plot_module, baseline_run in figures:
		dst = render_figure(name, run_key, plot_module, baseline_run, simouts)
		print('>>> figure {} -> {}'.format(
			name, os.path.relpath(dst, fp.ROOT_PATH)))

	if not args.no_metrics:
		print_metrics(simouts)


if __name__ == '__main__':
	main()
