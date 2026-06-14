"""Second-wave experiment — autocrine amplification of a weak transient agonist.

Under a saturating agonist the platelet Ca2+ response is store-limited, so the
v0.61 autocrine loops (secreted dense-granule ADP -> P2Y1; TXA2 -> TP -> Gq)
barely move the cytosolic Ca2+ trace (they show only in the receptor / messenger
readouts). Their effect on Ca2+ emerges at a WEAK, TRANSIENT stimulus: once the
primary agonist decays, PKC-triggered secretion releases dense-granule ADP that
sustains P2Y1, holding the store depleted so store-operated entry keeps cytosolic
Ca2+ elevated. This prolonged, autocrine-driven phase is the classic platelet
"second wave" -- and it is abolished in the open-loop (pre-v0.61) model.

This driver contrasts three closed/open-loop conditions at a weak ADP pulse
(ADP-only, +Ca2+, no thrombin / ATP) and writes the headline figure + data:

  v0.6      open loop    AUTOCRINE_ADP_GAIN=0, COX1_FACTOR=0  (pre-v0.61 baseline)
  aspirin   ADP loop on  AUTOCRINE_ADP_GAIN=1, COX1_FACTOR=0  (TXA2/COX-1 knocked out)
  full      both loops   AUTOCRINE_ADP_GAIN=1, COX1_FACTOR=1  (v0.61)

Decomposition: "full" vs "aspirin" isolates the TXA2 arm; "aspirin" vs "v0.6"
isolates the autocrine ADP arm. The sustained Ca2+ second wave is carried by the
ADP loop; TXA2 adds Gq drive that shows in IP3 but not Ca2+ (store/SOCE-limited).

Both knobs are restored after every run.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/runSecondWave.py [sim_outdir] \\
        [--adp-uM 0.5] [--length 300]

Outputs under out/<sim_outdir>/:
    second_wave_traces.png      3-panel (cyt Ca2+, IP3, P2Y1 desens) x 3 conditions
    second_wave.npz             time + per-condition trace matrices
    second_wave_summary.json    metadata + per-condition scalars
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from runscripts.manual.runPlateletSim import resolve_sim_path, run_platelet_sim
import reconstruction.platelet.dataclasses.process.calcium_signalling as cs_mod
import reconstruction.platelet.dataclasses.process.thromboxane_synthesis as tx_mod
import wholecell.utils.filepath as fp
from wholecell.io.tablereader import TableReader

@dataclass(frozen=True)
class Condition:
	"""Loop toggles + plot style/legend for one open/closed-loop condition."""
	adp_gain: float   # cs_mod.AUTOCRINE_ADP_GAIN
	cox1: float       # tx_mod.COX1_FACTOR
	label: str        # detailed legend entry
	color: str
	ls: str
	lw: float


CONDITIONS = {
	'v06': Condition(0.0, 0.0,
		'v0.6 — open loop (no autocrine feedback)', '#888888', '-', 2.6),
	'aspirin': Condition(1.0, 0.0,
		r'v0.61 + aspirin (TXA$_2$ / COX-1 knocked out; ADP loop only)',
		'#2980b9', '-.', 1.9),
	'full': Condition(1.0, 1.0,
		r'v0.61 — full model (autocrine ADP + TXA$_2$ loops)',
		'#c0392b', '--', 2.1),
}
ORDER = ['v06', 'aspirin', 'full']

# CalciumTrace columns harvested for the figure + scalars.
TRACE_COLS = ('ca_cyt_nM', 'ip3_nM', 'p2y1_desensitised_frac')


def run_condition(out_path: str, key: str, adp_uM: float, length: int,
		log_to_shell: bool = True) -> dict[str, np.ndarray]:
	"""Run one open/closed-loop condition; always restore the loop knobs."""
	cfg = CONDITIONS[key]
	adp0, cox0 = cs_mod.AUTOCRINE_ADP_GAIN, tx_mod.COX1_FACTOR
	cell_dir = os.path.join(out_path, f'{key}_cell')
	fp.makedirs(cell_dir)
	try:
		cs_mod.AUTOCRINE_ADP_GAIN = cfg.adp_gain
		tx_mod.COX1_FACTOR = cfg.cox1
		paths = run_platelet_sim(cell_dir, length_sec=length, seed=0,
			log_to_shell=False, ca_ex_mM=1.2,
			thrombin_peak_nM=0.0, adp_peak_uM=adp_uM, atp_ex_peak_uM=0.0)
	finally:
		cs_mod.AUTOCRINE_ADP_GAIN, tx_mod.COX1_FACTOR = adp0, cox0

	reader = TableReader(os.path.join(paths['sim_out_dir'], 'CalciumTrace'))
	traces = {c: reader.readColumn(c).flatten() for c in TRACE_COLS}
	if log_to_shell:
		cyt = traces['ca_cyt_nM']
		print(f'  {key:8s} peak={cyt.max():6.1f}  end={cyt[-1]:6.1f} nM  '
			f'(ADP gain {cfg.adp_gain:g}, COX-1 {cfg.cox1:g})')
	return traces


def reduce_scalars(traces: dict[str, np.ndarray]) -> dict[str, float]:
	cyt = traces['ca_cyt_nM']
	rest = float(cyt[0])
	t = np.arange(len(cyt), dtype=float)
	return {
		'rest_cyt_nM': rest,
		'peak_cyt_nM': float(cyt.max()),
		'end_cyt_nM': float(cyt[-1]),
		'sustained_auc_nMs': float(np.trapz(np.maximum(cyt - rest, 0.0), t)),
		'ip3_end_nM': float(traces['ip3_nM'][-1]),
		'p2y1_desens_max': float(traces['p2y1_desensitised_frac'].max()),
	}


def plot_second_wave(results: dict[str, dict[str, np.ndarray]],
		adp_uM: float, length: int, png_path: str) -> None:
	"""3-panel comparison with detailed legends + mathtext chemical formulae."""
	t = np.arange(length + 1)[: len(results['full']['ca_cyt_nM'])]
	panels = [
		('ca_cyt_nM', r'cytosolic $\mathrm{Ca^{2+}}$ (nM)',
			r'Sustained $\mathrm{Ca^{2+}}$ — the second wave'),
		('ip3_nM', r'$\mathrm{IP_3}$ (nM)',
			r'$\mathrm{IP_3}$ — autocrine $\mathrm{G_q}$ drive'),
		('p2y1_desensitised_frac', 'P2Y1 desensitised fraction',
			'P2Y1 engagement (PKC-phosphorylated)'),
	]
	fig, axes = plt.subplots(1, 3, figsize=(15, 4.9))
	for ax, (col, ylab, title) in zip(axes, panels):
		for key in ORDER:
			cfg = CONDITIONS[key]
			ax.plot(t, results[key][col], color=cfg.color,
				ls=cfg.ls, lw=cfg.lw, label=cfg.label)
		ax.set_xlabel('time (s)')
		ax.set_ylabel(ylab)
		ax.set_title(title, fontsize=10)
		ax.grid(alpha=0.3)
	axes[0].legend(frameon=False, fontsize=8, loc='upper right')

	# Annotate the sustained-Ca2+ gap (full vs open-loop) at the final time.
	cyt_full = results['full']['ca_cyt_nM']
	cyt_v06 = results['v06']['ca_cyt_nM']
	d_end = cyt_full[-1] - cyt_v06[-1]
	axes[0].annotate(
		f'+{d_end:.0f} nM\nsustained\n(autocrine ADP)',
		xy=(t[-1], cyt_full[-1]),
		xytext=(0.52 * t[-1], 0.45 * (cyt_full.max())),
		fontsize=8.5, color='#c0392b',
		arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.2))

	fig.suptitle(
		r'Weak transient agonist (ADP %g $\mu$M, ADP-only, +$\mathrm{Ca^{2+}}$ 1.2 mM, %d s): '
		'v0.61 autocrine loops sustain the response (second wave)'
		% (adp_uM, length), fontsize=12)
	# Detailed caption spelling out the mechanism + decomposition.
	fig.text(0.5, -0.02,
		'All conditions reach the same store-limited peak; they diverge in the '
		'recovery phase. Once the primary ADP pulse decays, the open-loop model '
		'(v0.6) returns toward baseline, while PKC-triggered dense-granule ADP '
		'secretion sustains P2Y1 in the closed-loop model, holding the store '
		'depleted so SOCE keeps cytosolic Ca2+ elevated. "full" vs "aspirin" '
		'isolates the TXA2 arm (visible in IP3, not Ca2+ — store/SOCE-limited); '
		'"aspirin" vs "v0.6" isolates the autocrine ADP arm, which carries the '
		'sustained Ca2+ second wave.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.02, 1, 0.95])
	fig.savefig(png_path, dpi=150, bbox_inches='tight')
	plt.close(fig)


def write_outputs(results, scalars, adp_uM, length, out_path):
	arrays = {'time': np.arange(length + 1)[: len(results['full']['ca_cyt_nM'])]}
	for key in ORDER:
		for col in TRACE_COLS:
			arrays[f'{key}__{col}'] = results[key][col]
	np.savez(os.path.join(out_path, 'second_wave.npz'), **arrays)

	plot_second_wave(results, adp_uM, length,
		os.path.join(out_path, 'second_wave_traces.png'))

	payload = {
		'description': 'Autocrine second-wave: open vs closed loop at a weak '
			'transient ADP stimulus.',
		'adp_uM': adp_uM, 'length_sec': length, 'ca_ex_mM': 1.2,
		'conditions': {key: {'toggles': {
			'AUTOCRINE_ADP_GAIN': CONDITIONS[key].adp_gain,
			'COX1_FACTOR': CONDITIONS[key].cox1},
			'scalars': scalars[key]} for key in ORDER},
		'sustained_cyt_gap_full_minus_v06_nM':
			scalars['full']['end_cyt_nM'] - scalars['v06']['end_cyt_nM'],
	}
	with open(os.path.join(out_path, 'second_wave_summary.json'), 'w') as f:
		json.dump(payload, f, indent='\t')


def run_second_wave(out_path: str, adp_uM: float = 0.5, length: int = 300,
		log_to_shell: bool = True):
	"""Run all three conditions, write figure + npz + summary; return results."""
	results, scalars = {}, {}
	for key in ORDER:
		results[key] = run_condition(out_path, key, adp_uM, length, log_to_shell)
		scalars[key] = reduce_scalars(results[key])
	write_outputs(results, scalars, adp_uM, length, out_path)
	return results, scalars


def _build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		description='Second-wave experiment: autocrine amplification of a weak '
			'transient agonist (open vs closed loop).')
	p.add_argument('sim_outdir', nargs='?', default=None,
		help='Output dir under out/. Default = second_wave_<timestamp>.')
	p.add_argument('--adp-uM', dest='adp_uM', type=float, default=0.5,
		help='Peak of the weak transient ADP pulse (µM). Default 0.5.')
	p.add_argument('--length', '--length-sec', dest='length_sec', type=int,
		default=300, help='Simulation length (s). Default 300.')
	p.add_argument('--no-log-to-shell', dest='log_to_shell',
		action='store_false', help='Suppress per-condition progress prints.')
	p.set_defaults(log_to_shell=True)
	return p


def main(argv: list[str] | None = None) -> None:
	args = _build_parser().parse_args(argv)
	sim_outdir = args.sim_outdir or 'second_wave_{}'.format(
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)
	print(f'Second-wave experiment (ADP {args.adp_uM} µM, {args.length_sec} s)')
	print(f'  Output: {sim_path}\n')
	_, scalars = run_second_wave(sim_path, adp_uM=args.adp_uM,
		length=args.length_sec, log_to_shell=args.log_to_shell)
	gap = scalars['full']['end_cyt_nM'] - scalars['v06']['end_cyt_nM']
	print(f'\n  Sustained Ca2+ gap (full - v0.6) at {args.length_sec} s: '
		f'{gap:+.1f} nM')
	print(f'  Wrote second_wave_traces.png / .npz / _summary.json to {sim_path}')


if __name__ == '__main__':
	main()
