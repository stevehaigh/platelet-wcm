"""
Phase 3 validation runner — Dolan 2014 Fig. 4 reproduction.

Runs the platelet calcium model under two conditions back-to-back:

  +Ca_ex   :  --ca-ex-mM 1.2  (Dolan nominal extracellular Ca²⁺)
  −Ca_ex   :  --ca-ex-mM 0    (Dolan EDTA condition)

Then invokes ``models.platelet.analysis.phase3_dolan_fig4.make_phase3_plot``
to produce a 3-panel comparison figure with PASS/FAIL annotations against
the Dolan filtering criteria. A JSON summary of the measured peaks /
plateaus / criteria results is written alongside the figure.

Output layout:

  out/<sim_outdir>/
      with_ca/         platelet_stub_<seed>/<seed>/generation_000000/000000/simOut/
      no_ca/           platelet_stub_<seed>/<seed>/generation_000000/000000/simOut/
      phase3_dolan_fig4.png
      phase3_summary.json
      kb/simData.cPickle             (shared)
      metadata/metadata.json         (run-level)

Usage:

  PYTHONPATH=$PWD python runscripts/manual/runPhase3.py [sim_outdir] \\
        [--length 200] [--seed 0] [--no-log-to-shell]

Default ``sim_outdir`` is ``phase3_<timestamp>``.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime

from models.platelet.analysis.phase3_dolan_fig4 import make_phase3_plot
from runscripts.manual.runPlateletSim import (
	resolve_sim_path,
	run_platelet_sim,
	write_metadata,
)
import wholecell.utils.filepath as fp


DEFAULT_LENGTH_SEC = 200      # Dolan Fig 4C uses 120-s window; 200 s gives margin


def _build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description='Run the Phase 3 Dolan 2014 Fig. 4 validation '
		            '(with vs without extracellular Ca²⁺).')
	parser.add_argument(
		'sim_outdir',
		nargs='?',
		default=None,
		help='Output directory under out/. Default = phase3_<timestamp>.')
	parser.add_argument(
		'--length', '--length-sec',
		dest='length_sec', type=int, default=DEFAULT_LENGTH_SEC,
		help='Simulation length in seconds for each condition. '
		     f'Default = {DEFAULT_LENGTH_SEC}.')
	parser.add_argument(
		'--seed', type=int, default=0,
		help='Random seed (shared across both conditions). Default = 0.')
	parser.add_argument(
		'--no-log-to-shell',
		dest='log_to_shell', action='store_false',
		help='Disable shell logging during the runs.')
	parser.set_defaults(log_to_shell=True)
	return parser


def _run_one_condition(condition_dir: str, label: str, length_sec: int,
		seed: int, ca_ex_mM: float, log_to_shell: bool) -> str:
	"""Run a single condition into condition_dir/, return its simOut path.

	Phase 3 always uses IP3 forcing (the spike is the stimulus we're
	testing the SOCE-shape response to); the comparison is over
	extracellular Ca²⁺ only.
	"""
	fp.makedirs(condition_dir)
	write_metadata(condition_dir,
		description=f'phase3 — {label}',
		seed=seed, length_sec=length_sec, ca_ex_mM=ca_ex_mM,
		ip3_forced=True)
	paths = run_platelet_sim(
		condition_dir,
		length_sec=length_sec,
		seed=seed,
		log_to_shell=log_to_shell,
		ca_ex_mM=ca_ex_mM,
		ip3_forced=True,
	)
	return paths['sim_out_dir']


def main(argv=None) -> None:
	args = _build_parser().parse_args(argv)

	sim_outdir = args.sim_outdir or 'phase3_{}'.format(
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)

	with_ca_dir = os.path.join(sim_path, 'with_ca')
	no_ca_dir   = os.path.join(sim_path, 'no_ca')

	print(f'Phase 3 — running both conditions into {sim_path}')
	print(f'  +Ca_ex (1.2 mM): {with_ca_dir}')
	print(f'  −Ca_ex (0):      {no_ca_dir}')
	print(f'  length: {args.length_sec} s, seed: {args.seed}')
	print()

	with_ca_simout = _run_one_condition(
		with_ca_dir, 'with_ca', args.length_sec, args.seed,
		ca_ex_mM=1.2, log_to_shell=args.log_to_shell)
	no_ca_simout = _run_one_condition(
		no_ca_dir, 'no_ca', args.length_sec, args.seed,
		ca_ex_mM=0.0, log_to_shell=args.log_to_shell)

	# ── Build the comparison figure and summary ────────────────────────────
	plot_path = os.path.join(sim_path, 'phase3_dolan_fig4.png')
	summary = make_phase3_plot(with_ca_simout, no_ca_simout, plot_path)

	summary_path = os.path.join(sim_path, 'phase3_summary.json')
	with open(summary_path, 'w') as f:
		json.dump(summary, f, indent='\t')

	# ── Console summary ────────────────────────────────────────────────────
	print()
	print('Phase 3 results:')
	print(f'  +Ca_ex: peak {summary["with_extracellular_ca"]["peak_cyt_nM"]:.1f} nM '
	      f'at t={summary["with_extracellular_ca"]["peak_t_s"]:.0f}s; '
	      f'DTS min {summary["with_extracellular_ca"]["dts_min_uM"]:.1f} µM')
	print(f'  −Ca_ex: peak {summary["without_extracellular_ca"]["peak_cyt_nM"]:.1f} nM '
	      f'at t={summary["without_extracellular_ca"]["peak_t_s"]:.0f}s; '
	      f'DTS min {summary["without_extracellular_ca"]["dts_min_uM"]:.1f} µM')
	print(f'  SOCE differential: {summary["soce_differential_nM"]:.1f} nM '
	      f'(criterion: ≥ 100 nM)')
	print()
	passed = sum(1 for c in summary['criteria'] if c['passed'])
	total  = len(summary['criteria'])
	print(f'Acceptance criteria: {passed}/{total} pass')
	for c in summary['criteria']:
		mark = '✓' if c['passed'] else '✗'
		print(f'  {mark} {c["criterion"]:<30}  {c["rule"]:<32}  → {c["measured"]}')
	print()
	print(f'Figure: {plot_path}')
	print(f'Summary: {summary_path}')


if __name__ == '__main__':
	main()
