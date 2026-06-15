"""Integrin αIIbβ3 inside-out activation figure (v0.61 §3).

PKC + CalDAG-GEFI/Ca²⁺ → Rap1b → talin switches αIIbβ3 (GPIIb-IIIa) from its
resting (low-affinity) to its active (high-affinity) conformation. The active
fraction is the per-cell PAC-1 flow-cytometry readout. This figure shows:

  * graded PAC-1 activation across agonist conditions (standard > thrombin-only),
  * exactly zero at rest and under the αIIbβ3-antagonist / Glanzmann knockout
    (integrin_act_scale = 0) — even though the PKC × Ca²⁺ gate still fires,
  * the resting ⇌ active conformational switch (mass-conserving) under the
    standard transient.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotIntegrin.py [sim_outdir] \\
        [--length 250]

Outputs under out/<sim_outdir>/:
    integrin_activation.png   3-panel (PAC-1 fraction, gate, state populations)
    integrin_activation.npz   time + per-condition trace matrices
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field
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
	"""RunConfig knobs + plot style for one integrin condition."""
	label: str
	color: str
	ls: str
	kwargs: dict = field(default_factory=dict)


# All conditions are +Ca²⁺ (1.2 mM). kwargs override agonist peaks / the
# integrin knockout knob relative to the standard transient.
CONDITIONS = {
	'standard': Condition(
		'standard transient (thrombin + ADP + ATP)', '#c0392b', '-', {}),
	'thrombin': Condition(
		'thrombin-only', '#e67e22', '-',
		dict(adp_peak_uM=0.0, atp_ex_peak_uM=0.0)),
	'adp': Condition(
		'ADP-only', '#8e44ad', '-.',
		dict(thrombin_peak_nM=0.0, atp_ex_peak_uM=0.0)),
	'glanzmann': Condition(
		r'Glanzmann / αIIbβ3 antagonist (scale = 0)', '#2980b9', '--',
		dict(integrin_act_scale=0.0)),
	'rest': Condition(
		'resting (no agonist)', '#7f8c8d', ':',
		dict(thrombin_peak_nM=0.0, adp_peak_uM=0.0, atp_ex_peak_uM=0.0)),
}
ORDER = ['standard', 'thrombin', 'adp', 'glanzmann', 'rest']

TRACE_COLS = ('active_frac', 'integrin_gate', 'aIIbb3_active', 'aIIbb3_resting')


def run_condition(out_path: str, key: str, length: int,
		log_to_shell: bool = True) -> dict[str, np.ndarray]:
	cfg = CONDITIONS[key]
	cell_dir = os.path.join(out_path, f'{key}_cell')
	fp.makedirs(cell_dir)
	run_config = RunConfig(ca_ex_mM=1.2, **cfg.kwargs)
	paths = run_platelet_sim(cell_dir, length_sec=length, seed=0,
		log_to_shell=False, run_config=run_config)
	reader = TableReader(os.path.join(paths['sim_out_dir'], 'IntegrinTrace'))
	traces = {c: reader.readColumn(c).flatten() for c in TRACE_COLS}
	if log_to_shell:
		print(f'  {key:10s} PAC-1 active_frac: peak={traces["active_frac"].max():.3f} '
			f'end={traces["active_frac"][-1]:.3f}')
	return traces


def plot_integrin(results: dict[str, dict[str, np.ndarray]], length: int,
		png_path: str) -> None:
	t = np.arange(length + 1)[: len(results['standard']['active_frac'])]
	fig, axes = plt.subplots(1, 3, figsize=(15, 4.9))

	# Panel 1 — PAC-1 active fraction across conditions.
	for key in ORDER:
		cfg = CONDITIONS[key]
		axes[0].plot(t, results[key]['active_frac'], color=cfg.color,
			ls=cfg.ls, lw=2.0, label=cfg.label)
	axes[0].set_ylabel(r'$\alpha_{IIb}\beta_3$ active fraction (PAC-1$^+$)')
	axes[0].set_title('PAC-1 — high-affinity integrin', fontsize=10)
	axes[0].legend(frameon=False, fontsize=7.5, loc='upper left')

	# Panel 2 — the inside-out activation gate (PKC × Ca²⁺).
	for key in ORDER:
		cfg = CONDITIONS[key]
		axes[1].plot(t, results[key]['integrin_gate'], color=cfg.color,
			ls=cfg.ls, lw=2.0, label=cfg.label)
	axes[1].set_ylabel(r'inside-out gate (PKC $\times$ $\mathrm{Ca^{2+}}$)')
	axes[1].set_title('Activation gate (the driver)', fontsize=10)
	axes[1].annotate('Glanzmann: gate fires,\nbut no activation',
		xy=(0.5, 0.5), xycoords='axes fraction', fontsize=8,
		color='#2980b9', ha='center')

	# Panel 3 — the resting ⇌ active conformational switch (standard).
	std = results['standard']
	total = std['aIIbb3_active'] + std['aIIbb3_resting']
	axes[2].plot(t, std['aIIbb3_resting'] / 1000.0, color='#16a085', lw=2.0,
		label='resting (low-affinity)')
	axes[2].plot(t, std['aIIbb3_active'] / 1000.0, color='#c0392b', lw=2.0,
		label='active (high-affinity)')
	axes[2].plot(t, total / 1000.0, color='#7f8c8d', ls=':', lw=1.5,
		label='total (conserved)')
	axes[2].set_ylabel(r'$\alpha_{IIb}\beta_3$ count ($\times 10^3$)')
	axes[2].set_title('Conformational switch (standard transient)', fontsize=10)
	axes[2].legend(frameon=False, fontsize=7.5, loc='center right')

	for ax in axes:
		ax.set_xlabel('time (s)')
		ax.grid(alpha=0.3)

	fig.suptitle(
		r'Integrin $\alpha_{IIb}\beta_3$ inside-out activation (v0.61 §3) — '
		'the terminal PKC output', fontsize=12)
	fig.text(0.5, -0.04,
		'PKC + CalDAG-GEFI/Ca²⁺ → Rap1b → talin switches αIIbβ3 (GPIIb-IIIa) from '
		'low- to high-affinity; the active fraction is the per-cell PAC-1 readout '
		'(flow cytometry). Left: PAC-1 activation is graded by agonist strength '
		'(standard > thrombin-only > ADP-only), reaches ~0.77 under the standard '
		'transient, and is exactly zero at rest and under the αIIbβ3-antagonist / '
		'Glanzmann knockout. Middle: the shared PKC × Ca²⁺ coincidence gate drives '
		'activation — note it still fires in the Glanzmann case (the receptor is '
		'blocked downstream, not the signal), and is exactly zero at rest. Right: '
		'the resting ⇌ active conformational switch is mass-conserving (total '
		'αIIbβ3 constant). Aggregation itself is inter-cellular and out of scope '
		'for a single-platelet model — only the affinity state is represented.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.04, 1, 0.95])
	fig.savefig(png_path, dpi=150, bbox_inches='tight')
	plt.close(fig)


def run_integrin_figure(out_path: str, length: int = 250,
		log_to_shell: bool = True):
	results = {}
	for key in ORDER:
		results[key] = run_condition(out_path, key, length, log_to_shell)
	arrays = {'time': np.arange(length + 1)[: len(results['standard']['active_frac'])]}
	for key in ORDER:
		for col in TRACE_COLS:
			arrays[f'{key}__{col}'] = results[key][col]
	np.savez(os.path.join(out_path, 'integrin_activation.npz'), **arrays)
	plot_integrin(results, length,
		os.path.join(out_path, 'integrin_activation.png'))
	return results


def _build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		description='Integrin αIIbβ3 inside-out activation figure (v0.61 §3).')
	p.add_argument('sim_outdir', nargs='?', default=None,
		help='Output dir under out/. Default = integrin_<timestamp>.')
	p.add_argument('--length', '--length-sec', dest='length_sec', type=int,
		default=250, help='Simulation length (s). Default 250.')
	return p


def main(argv: list[str] | None = None) -> None:
	args = _build_parser().parse_args(argv)
	sim_outdir = args.sim_outdir or 'integrin_{}'.format(
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)
	print(f'Integrin αIIbβ3 activation figure ({args.length_sec} s)')
	print(f'  Output: {sim_path}\n')
	run_integrin_figure(sim_path, length=args.length_sec)
	print(f'\n  Wrote integrin_activation.png / .npz to {sim_path}')


if __name__ == '__main__':
	main()
