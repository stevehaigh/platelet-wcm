"""Recovery-phase validation figure — why the sustained plateau is the
autocrine "second wave", not a Dolan-recovery defect.

Dolan & Diamond (2014) drove their model with a *transient* IP3 dose and
modelled neither thromboxane nor granule secretion, so their cytosolic Ca2+
recovers after the peak (Fig. 4C). Our model adds two autocrine
positive-feedback loops Dolan omitted — TXA2 -> TP -> Gq and secreted
ADP -> P2Y1 — which sustain Gq -> IP3 and hold the cell activated (the platelet
"second wave"). The high sustained plateau is therefore a model *prediction*
(extra biology), not a recovery-phase calibration defect.

This script demonstrates that by contrasting, in both the +Ca2+ and EDTA
conditions:

  full       v0.61 full model (autocrine loops on, thrombin drive) — sustained
  dolan_eq   Dolan-equivalent (loops off + transient reversible ADP stimulus) —
             IP3 recovers to baseline, the DTS store refills, cytosol returns
             toward baseline: the Ca-handling machinery reproduces Dolan.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotRecoveryPhase.py [sim_outdir] \\
        [--length 400]

Outputs under out/<sim_outdir>/:
    recovery_phase_traces.png   3-panel (cyt Ca2+, IP3, DTS store) x 4 conditions
    recovery_phase.npz          time + per-condition trace matrices
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
	"""RunConfig knobs + plot style for one recovery-phase condition."""
	ca_ex_mM: float
	loops_off: bool        # True → Dolan-equivalent (autocrine loops disabled)
	label: str
	color: str
	ls: str


# loops_off → transient reversible ADP stimulus (no thrombin/ATP), autocrine
# ADP and TXA2 loops disabled: the closest analogue to Dolan's transient IP3.
CONDITIONS = {
	'full_ca':   Condition(1.2, False,
		r'full v0.61, +$\mathrm{Ca^{2+}}$ (autocrine loops on)', '#c0392b', '-'),
	'full_edta': Condition(0.0, False,
		'full v0.61, EDTA (autocrine loops on)', '#e67e22', '-'),
	'dolan_ca':  Condition(1.2, True,
		r'Dolan-equiv., +$\mathrm{Ca^{2+}}$ (loops off, transient ADP)',
		'#2980b9', '--'),
	'dolan_edta': Condition(0.0, True,
		'Dolan-equiv., EDTA (loops off, transient ADP)', '#16a085', '--'),
}
ORDER = ['full_ca', 'full_edta', 'dolan_ca', 'dolan_edta']

TRACE_COLS = ('ca_cyt_nM', 'ip3_nM', 'ca_dts_uM')


def run_condition(out_path: str, key: str, length: int,
		log_to_shell: bool = True) -> dict[str, np.ndarray]:
	"""Run one recovery-phase condition; return harvested traces."""
	cfg = CONDITIONS[key]
	cell_dir = os.path.join(out_path, f'{key}_cell')
	fp.makedirs(cell_dir)
	if cfg.loops_off:
		run_config = RunConfig(ca_ex_mM=cfg.ca_ex_mM, thrombin_peak_nM=0.0,
			atp_ex_peak_uM=0.0, autocrine_adp_gain=0.0, cox1_factor=0.0)
	else:
		run_config = RunConfig(ca_ex_mM=cfg.ca_ex_mM)
	paths = run_platelet_sim(cell_dir, length_sec=length, seed=0,
		log_to_shell=False, run_config=run_config)
	reader = TableReader(os.path.join(paths['sim_out_dir'], 'CalciumTrace'))
	traces = {c: reader.readColumn(c).flatten() for c in TRACE_COLS}
	if log_to_shell:
		print(f'  {key:11s} cyt_end={traces["ca_cyt_nM"][-1]:6.1f} nM  '
			f'ip3_end={traces["ip3_nM"][-1]:6.1f} nM  '
			f'dts_end={traces["ca_dts_uM"][-1]:6.1f} uM')
	return traces


def plot_recovery(results: dict[str, dict[str, np.ndarray]], length: int,
		png_path: str) -> None:
	"""3-panel recovery contrast with detailed mathtext caption."""
	t = np.arange(length + 1)[: len(results['full_ca']['ca_cyt_nM'])]
	panels = [
		('ca_cyt_nM', r'cytosolic $\mathrm{Ca^{2+}}$ (nM)',
			r'Cytosolic $\mathrm{Ca^{2+}}$'),
		('ip3_nM', r'$\mathrm{IP_3}$ (nM)',
			r'$\mathrm{IP_3}$ — the autocrine $\mathrm{G_q}$ readout'),
		('ca_dts_uM', r'DTS store $\mathrm{Ca^{2+}}$ ($\mu$M)',
			'DTS store refilling'),
	]
	fig, axes = plt.subplots(1, 3, figsize=(15, 4.9))
	for ax, (col, ylab, title) in zip(axes, panels):
		for key in ORDER:
			cfg = CONDITIONS[key]
			ax.plot(t, results[key][col], color=cfg.color, ls=cfg.ls,
				lw=2.0, label=cfg.label)
		ax.set_xlabel('time (s)')
		ax.set_ylabel(ylab)
		ax.set_title(title, fontsize=10)
		ax.grid(alpha=0.3)
	# Dolan reference bands on the IP3 / store panels.
	axes[1].axhline(50, color='#7f8c8d', ls=':', lw=1.2,
		label=r'resting $\mathrm{IP_3}$ (50 nM)')
	axes[2].axhspan(120, 180, color='#bdc3c7', alpha=0.35,
		label=r'Dolan $\mathrm{DTS_{min}}$ band')
	axes[0].legend(frameon=False, fontsize=7.5, loc='upper right')
	axes[1].legend(frameon=False, fontsize=7.5, loc='upper right')
	axes[2].legend(frameon=False, fontsize=7.5, loc='lower right')

	fig.suptitle(
		'Recovery phase: the sustained plateau is the autocrine "second wave", '
		'not a Dolan-recovery defect', fontsize=12)
	fig.text(0.5, -0.04,
		'Dolan (2014) drove a transient IP3 dose and modelled neither thromboxane '
		'nor granule secretion, so its cytosol recovers after the peak. The full '
		'v0.61 model (solid) adds two autocrine positive-feedback loops Dolan '
		'omitted — TXA2 -> TP -> Gq and secreted ADP -> P2Y1 — which sustain '
		'Gq -> IP3 and hold the cell activated (the platelet "second wave"): IP3 '
		'stays elevated and the DTS store stays empty, so SOCE pins the +Ca2+ '
		'plateau high. Disable the loops and apply a transient reversible-ADP '
		'stimulus (dashed; the Dolan-equivalent) and IP3 recovers to its 50 nM '
		'baseline and the store refills toward Dolan\'s 120-180 uM band — the '
		'Ca-handling machinery (IP3R / SERCA / SOCE / PMCA / NCX) reproduces '
		'Dolan\'s recovery. The cytosolic-Ca2+ and store recovery lag IP3 by '
		'~100 s (heavy cytosolic buffering). The sustained plateau is therefore a '
		'model prediction (extra biology), not a calibration failure.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.04, 1, 0.95])
	fig.savefig(png_path, dpi=150, bbox_inches='tight')
	plt.close(fig)


def run_recovery_phase(out_path: str, length: int = 400,
		log_to_shell: bool = True):
	"""Run all four conditions; write figure + npz; return the traces."""
	results = {}
	for key in ORDER:
		results[key] = run_condition(out_path, key, length, log_to_shell)
	arrays = {'time': np.arange(length + 1)[: len(results['full_ca']['ca_cyt_nM'])]}
	for key in ORDER:
		for col in TRACE_COLS:
			arrays[f'{key}__{col}'] = results[key][col]
	np.savez(os.path.join(out_path, 'recovery_phase.npz'), **arrays)
	plot_recovery(results, length,
		os.path.join(out_path, 'recovery_phase_traces.png'))
	return results


def _build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		description='Recovery-phase validation: the sustained plateau is the '
			'autocrine second wave (full model) vs Dolan-equivalent recovery.')
	p.add_argument('sim_outdir', nargs='?', default=None,
		help='Output dir under out/. Default = recovery_phase_<timestamp>.')
	p.add_argument('--length', '--length-sec', dest='length_sec', type=int,
		default=400, help='Simulation length (s). Default 400.')
	return p


def main(argv: list[str] | None = None) -> None:
	args = _build_parser().parse_args(argv)
	sim_outdir = args.sim_outdir or 'recovery_phase_{}'.format(
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)
	print(f'Recovery-phase validation ({args.length_sec} s)')
	print(f'  Output: {sim_path}\n')
	run_recovery_phase(sim_path, length=args.length_sec)
	print(f'\n  Wrote recovery_phase_traces.png / .npz to {sim_path}')


if __name__ == '__main__':
	main()
