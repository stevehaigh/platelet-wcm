"""MCU-knockout figure — with the #76 coupling, MCU loss *reduces* the
agonist-evoked cytosolic Ca2+ peak (Ghatge 2026 / Ajanel 2025 direction), and
toggling the coupling off recovers the superseded buffer-only divergence.

This regenerates the thesis Results MCU figure for the *coupled* v0.6 model,
replacing the v0.5 buffer-only snapshot. Three conditions under a standard +Ca2+
agonist transient (60 s settle, 300 s), all else equal:

  wt        wild type, MCU uptake coupled to IP3R relief (mito_coupling_gain=1)
  ko        MCU knockout (mcu_vmax_scale=0), coupling on  — peak REDUCED
  ko_decpl  MCU knockout, coupling OFF (mito_coupling_gain=0) — buffer-only,
            peak RAISED (the superseded behaviour the coupling corrects)

The coupling (ip3r_relief_factor): MCU uptake at mitochondria-DTS (MAM) contact
sites relieves the IP3R's Ca2+-dependent inactivation, so losing MCU reduces
agonist-evoked release rather than merely redistributing cytosolic Ca2+. See
reports/design/mcu-coupling-2026-06-23.qmd.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotMcuCoupling.py [sim_outdir] \\
        [--length 300]

Outputs under out/<sim_outdir>/:
    mcu_coupling_traces.png   2-panel (cyt Ca2+, DTS store) x 3 conditions
    mcu_coupling.npz          time + per-condition trace matrices
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from runscripts.manual.runPlateletSim import resolve_sim_path, run_platelet_sim
from reconstruction.platelet.run_config import RunConfig
import wholecell.utils.filepath as fp
from wholecell.io.tablereader import TableReader


@dataclass(frozen=True)
class Condition:
	"""RunConfig knobs + plot style for one MCU condition."""
	mcu_vmax_scale: float
	mito_coupling_gain: float
	label: str
	color: str
	ls: str


CONDITIONS = {
	'wt': Condition(1.0, 1.0,
		'wild type (MCU coupled)', '#222222', '-'),
	'ko': Condition(0.0, 1.0,
		r'MCU knockout (coupled) — peak $\downarrow$', '#c0392b', '-'),
	'ko_decpl': Condition(0.0, 0.0,
		r'MCU knockout (decoupled, buffer-only) — peak $\uparrow$',
		'#7f8c8d', '--'),
}
ORDER = ['wt', 'ko', 'ko_decpl']

TRACE_COLS = ('ca_cyt_nM', 'ca_dts_uM')

# Standard +Ca2+ agonist transient with a 60 s settle — matches the acceptance /
# validation-target MCU runs (WT ~530, KO ~435, decoupled KO ~602 nM).
CA_EX_MM = 1.2
AGONIST_DELAY_S = 60.0


def run_condition(out_path: str, key: str, length: int,
		log_to_shell: bool = True) -> dict[str, np.ndarray]:
	"""Run one MCU condition; return harvested CalciumTrace columns."""
	cfg = CONDITIONS[key]
	cell_dir = os.path.join(out_path, f'{key}_cell')
	fp.makedirs(cell_dir)
	run_config = RunConfig(ca_ex_mM=CA_EX_MM, agonist_delay_s=AGONIST_DELAY_S,
		mcu_vmax_scale=cfg.mcu_vmax_scale,
		mito_coupling_gain=cfg.mito_coupling_gain)
	paths = run_platelet_sim(cell_dir, length_sec=length, seed=0,
		log_to_shell=False, run_config=run_config)
	reader = TableReader(os.path.join(paths['sim_out_dir'], 'CalciumTrace'))
	traces = {c: reader.readColumn(c).flatten() for c in TRACE_COLS}
	if log_to_shell:
		print(f'  {key:9s} peak_cyt={traces["ca_cyt_nM"].max():6.1f} nM  '
			f'dts_min={traces["ca_dts_uM"].min():6.2f} uM')
	return traces


def plot_mcu_coupling(results: dict[str, dict[str, np.ndarray]], length: int,
		png_path: str) -> None:
	"""2-panel cyt + DTS contrast with a detailed mathtext caption."""
	t = np.arange(length + 1)[: len(results['wt']['ca_cyt_nM'])]
	wt_peak = float(results['wt']['ca_cyt_nM'].max())
	ko_peak = float(results['ko']['ca_cyt_nM'].max())
	dec_peak = float(results['ko_decpl']['ca_cyt_nM'].max())

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
	for key in ORDER:
		cfg = CONDITIONS[key]
		ax1.plot(t, results[key]['ca_cyt_nM'], color=cfg.color, ls=cfg.ls,
			lw=2.0, label=cfg.label)
		ax2.plot(t, results[key]['ca_dts_uM'], color=cfg.color, ls=cfg.ls,
			lw=2.0, label=cfg.label)
	ax1.set_xlabel('time (s)')
	ax1.set_ylabel(r'cytosolic $\mathrm{Ca^{2+}}$ (nM)')
	ax1.set_title(r'MCU loss reduces the evoked cytosolic $\mathrm{Ca^{2+}}$ peak')
	ax1.legend(frameon=False, fontsize=8.5, loc='upper left')
	ax1.grid(alpha=0.3)

	ax2.set_xlabel('time (s)')
	ax2.set_ylabel(r'DTS store $\mathrm{Ca^{2+}}$ ($\mu$M)')
	ax2.set_title('the DTS store depletes (least for the coupled knockout)')
	ax2.legend(frameon=False, fontsize=8.5, loc='upper right')
	ax2.grid(alpha=0.3)

	fig.suptitle('MCU knockout and the mitochondria–DTS coupling '
		r'(+$\mathrm{Ca^{2+}}$, ' + f'{length:.0f} s)', fontsize=12)
	fig.text(0.5, -0.05,
		'With mitochondrial uptake coupled to IP3R release at mitochondria-DTS '
		'(MAM) contacts (ip3r_relief_factor; the shipped default), MCU knockout '
		f'(red) LOWERS the evoked cytosolic peak ({wt_peak:.0f} -> {ko_peak:.0f} '
		'nM), reproducing the direction of platelet MCU-knockout data (Ghatge '
		'2026; Ajanel 2025). Disabling the coupling (grey dashed; '
		'mito_coupling_gain=0) recovers the superseded buffer-only behaviour, in '
		f'which the same knockout instead RAISES the peak ({dec_peak:.0f} nM) — '
		'the divergence the coupling was added to correct. The effect is modest '
		'because losing the mitochondrial buffer partly offsets the lost relief. '
		'The DTS store (right) depletes in every case, least for the coupled '
		'knockout — the reduced IP3R release that drives its lower cytosolic '
		'peak. Deterministic model, one run per condition; seed 0.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.05, 1, 0.95])
	fig.savefig(png_path, dpi=150, bbox_inches='tight')
	plt.close(fig)


def run_mcu_coupling(out_path: str, length: int = 300,
		log_to_shell: bool = True) -> dict[str, dict[str, np.ndarray]]:
	"""Run all three conditions; write figure + npz; return the traces."""
	results = {}
	for key in ORDER:
		results[key] = run_condition(out_path, key, length, log_to_shell)
	arrays = {'time': np.arange(length + 1)[: len(results['wt']['ca_cyt_nM'])]}
	for key in ORDER:
		for col in TRACE_COLS:
			arrays[f'{key}__{col}'] = results[key][col]
	np.savez(os.path.join(out_path, 'mcu_coupling.npz'), **arrays)
	plot_mcu_coupling(results, length,
		os.path.join(out_path, 'mcu_coupling_traces.png'))
	return results


def _build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		description='MCU-knockout figure: with the #76 coupling, MCU loss '
			'reduces the evoked cytosolic Ca2+ peak; decoupling recovers the '
			'buffer-only divergence.')
	p.add_argument('sim_outdir', nargs='?', default=None,
		help='Output dir under out/. Default = mcu_coupling_<timestamp>.')
	p.add_argument('--length', '--length-sec', dest='length_sec', type=int,
		default=300, help='Simulation length (s). Default 300.')
	return p


def main(argv: list[str] | None = None) -> None:
	args = _build_parser().parse_args(argv)
	sim_outdir = args.sim_outdir or 'mcu_coupling_{}'.format(
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)
	print(f'MCU-coupling figure ({args.length_sec} s)')
	print(f'  Output: {sim_path}\n')
	run_mcu_coupling(sim_path, length=args.length_sec)
	print(f'\n  Wrote mcu_coupling_traces.png / .npz to {sim_path}')


if __name__ == '__main__':
	main()
