"""Single-mechanism perturbation experiments (PMCA V_max, MCU uptake).

Two in-silico knockdown/knockout scans on the frozen calcium model, producing
the mechanistic-finding figures the thesis Results section claims (issue #53).
Each condition overrides one rate constant in `calcium_signalling`, runs a
platelet sim, harvests the cytosolic / DTS Ca²⁺ traces, then restores the
constant — the same live-override pattern `runDoseSweep.py` uses.

Exp A — PMCA V_max rate-limits cytosolic Ca²⁺ recovery
    Knob   : K_PMCA['k_cat']  (basal PMCA turnover, baseline 5.5 s⁻¹)
    Scan   : ×{0.25, 0.5, 1, 2, 4}
    Protocol: EDTA / no-extracellular-Ca²⁺ (ca_ex = 0). SOCE cannot refill the
        store, so the DTS releases one bolus and PMCA/NCX clear the cytosol — a
        self-limiting transient whose recovery rate PMCA sets. (Under +Ca²⁺ the
        thrombin-cleaved PARs keep the cell activated and cyt never returns to
        rest, so recovery is undefined — see the lab note for issue #53.)
    Observable: recovery-tail AUC = ∫ max(cyt − rest, 0) dt over the post-peak
        window. Censoring-free; decreases as PMCA V_max rises.

Exp B — MCU buffers cytosolic Ca²⁺ without rescuing the DTS store
    Knob   : K_MITO['V_max_MCU']  (baseline 50000 ions/s; 0 = knockout)
    Scan   : {0 (KO), 1 (baseline), 4 (over-expression)}
    Protocol: +Ca²⁺ (ca_ex = 1.2 mM), standard sustained agonist transient.
    Observable: peak/sustained cyt (MCU lowers it — buffering) and DTS depletion
        (full in every case — MCU doesn't refill the store). MCU=0 reproduces
        the elevated-cytosolic-Ca²⁺ phenotype of MCU⁻/⁻ platelets (Ghatge 2026).

Usage:
    PYTHONPATH=$PWD python runscripts/manual/runPerturbation.py [sim_outdir] \\
        [--experiment pmca|mcu|both] [--length N] [--keep-cell-output]

Outputs under out/<sim_outdir>/:
    <exp>_traces.png             overlaid traces, colour-graded by factor
    <exp>.npz                    factors + cyt/DTS trace matrices + scalars
    perturbation_summary.json    metadata + per-condition scalars
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from dataclasses import dataclass, field
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from runscripts.manual.runPlateletSim import resolve_sim_path, run_platelet_sim
from reconstruction.platelet.dataclasses.process import calcium_signalling as cs_mod
import wholecell.utils.filepath as fp
from wholecell.io.tablereader import TableReader


# ── Experiment registry ─────────────────────────────────────────────────────
# Each experiment overrides ONE constant: cs_mod.<dict_name>[<knob_key>] is
# multiplied by each factor in turn. ca_ex_mM and length_sec set the protocol.
EXPERIMENTS: dict[str, dict] = {
	'pmca': dict(
		title='PMCA V_max rate-limits cytosolic Ca recovery (EDTA transient)',
		dict_name='K_PMCA', knob_key='k_cat', knob_label='PMCA V_max',
		factors=[0.25, 0.5, 1.0, 2.0, 4.0],
		ca_ex_mM=0.0, length_sec=300,
	),
	'mcu': dict(
		title='MCU buffers cytosolic Ca without rescuing the DTS store (+Ca)',
		dict_name='K_MITO', knob_key='V_max_MCU', knob_label='MCU V_max',
		factors=[0.0, 1.0, 4.0],
		ca_ex_mM=1.2, length_sec=400,
	),
}


def harvest_traces(simout_dir: str) -> tuple[np.ndarray, np.ndarray]:
	"""Read the cytosolic (nM) and DTS (µM) Ca²⁺ time series from CalciumTrace."""
	reader = TableReader(os.path.join(simout_dir, 'CalciumTrace'))
	cyt = reader.readColumn('ca_cyt_nM').flatten()
	dts = reader.readColumn('ca_dts_uM').flatten()
	return cyt, dts


def reduce_scalars(cyt: np.ndarray, dts: np.ndarray) -> dict[str, float]:
	"""Reduce a condition's traces to comparison scalars."""
	t = np.arange(len(cyt), dtype=float)        # 1-s timesteps
	rest = float(cyt[0])
	peak_idx = int(cyt.argmax())
	tail = np.maximum(cyt[peak_idx:] - rest, 0.0)
	return {
		'rest_cyt_nM':          rest,
		'peak_cyt_nM':          float(cyt.max()),
		't_peak_s':             float(peak_idx),
		'cyt_end_nM':           float(cyt[-1]),
		'recovery_tail_auc_nMs': float(np.trapz(tail, t[peak_idx:])),
		'dts_rest_uM':          float(dts[0]),
		'dts_min_uM':           float(dts.min()),
		'dts_end_uM':           float(dts[-1]),
	}


@dataclass
class PerturbationScan:
	"""Outcome of a single-knob perturbation scan."""
	exp_key: str
	cfg: dict
	baseline_value: float
	length_sec: int
	factors: list[float]
	cyt: np.ndarray                              # (n_factor, n_time)
	dts: np.ndarray                              # (n_factor, n_time)
	scalars: list[dict] = field(default_factory=list)

	def to_npz(self, path: str) -> None:
		np.savez(path,
			factors=np.array(self.factors), cyt=self.cyt, dts=self.dts,
			baseline_value=self.baseline_value, length_sec=self.length_sec)


# ── Driver ────────────────────────────────────────────────────────────────

def run_perturbation(out_path: str, exp_key: str,
		length_override: int | None = None,
		keep_cell_output: bool = False,
		log_to_shell: bool = True) -> PerturbationScan:
	"""Run one experiment's factor scan; always restore the knob afterwards."""
	cfg = EXPERIMENTS[exp_key]
	knob = getattr(cs_mod, cfg['dict_name'])
	key = cfg['knob_key']
	baseline = float(knob[key])
	length = length_override or cfg['length_sec']
	factors = list(cfg['factors'])

	cells_root = os.path.join(out_path, f'{exp_key}_cells')
	fp.makedirs(cells_root)

	cyt_rows: list[np.ndarray] = []
	dts_rows: list[np.ndarray] = []
	scalars: list[dict] = []
	try:
		for f in factors:
			knob[key] = baseline * f
			cell_dir = os.path.join(cells_root, f'x{f:g}')
			fp.makedirs(cell_dir)
			paths = run_platelet_sim(cell_dir, length_sec=length, seed=0,
				log_to_shell=False, ca_ex_mM=cfg['ca_ex_mM'])
			cyt, dts = harvest_traces(paths['sim_out_dir'])
			cyt_rows.append(cyt)
			dts_rows.append(dts)
			sc = {'factor': float(f), **reduce_scalars(cyt, dts)}
			scalars.append(sc)
			if log_to_shell:
				print(f'  {cfg["knob_label"]} ×{f:<5g} → '
					f'peak cyt {sc["peak_cyt_nM"]:>6.1f} nM, '
					f'recovery-AUC {sc["recovery_tail_auc_nMs"]:>9.0f} nM·s, '
					f'DTS min {sc["dts_min_uM"]:>6.1f} µM')
			if not keep_cell_output:
				_prune_cell(cell_dir)
	finally:
		knob[key] = baseline   # restore no matter what

	return PerturbationScan(
		exp_key=exp_key, cfg=cfg, baseline_value=baseline, length_sec=length,
		factors=factors, cyt=np.array(cyt_rows), dts=np.array(dts_rows),
		scalars=scalars)


def _prune_cell(cell_dir: str) -> None:
	"""Keep only simOut/CalciumTrace/ for each condition (mirrors dose sweep)."""
	for root, _, _ in os.walk(cell_dir):
		if os.path.basename(root) == 'simOut':
			for entry in os.listdir(root):
				if entry == 'CalciumTrace':
					continue
				full = os.path.join(root, entry)
				shutil.rmtree(full) if os.path.isdir(full) else os.remove(full)
	for d in ('kb', 'metadata'):
		full = os.path.join(cell_dir, d)
		if os.path.isdir(full):
			shutil.rmtree(full)


# ── Plots ─────────────────────────────────────────────────────────────────

def plot_pmca(scan: PerturbationScan, png_path: str) -> None:
	"""2-panel: cyt recovery traces + recovery-tail AUC vs PMCA V_max."""
	t = np.arange(scan.cyt.shape[1])
	cmap = plt.get_cmap('viridis')
	norm = np.linspace(0, 1, len(scan.factors))

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
	for i, f in enumerate(scan.factors):
		ax1.plot(t, scan.cyt[i], color=cmap(norm[i]), linewidth=1.8,
			label=f'×{f:g}')
	ax1.set_xlabel('Time (s)')
	ax1.set_ylabel(r'Cytosolic Ca$^{2+}$ (nM)')
	ax1.set_title('Cytosolic recovery vs PMCA V$_{max}$ (EDTA)')
	ax1.legend(title='PMCA V$_{max}$', fontsize=9)
	ax1.grid(alpha=0.3)

	auc = [s['recovery_tail_auc_nMs'] for s in scan.scalars]
	ax2.plot(scan.factors, auc, 'o-', color='#0b5394', linewidth=1.8, markersize=7)
	ax2.set_xscale('log')
	ax2.set_xlabel(r'PMCA V$_{max}$ (× baseline, log)')
	ax2.set_ylabel(r'Recovery-tail AUC (nM$\cdot$s)')
	ax2.set_title('Integrated residual Ca$^{2+}$ vs PMCA V$_{max}$')
	ax2.grid(alpha=0.3, which='both')

	fig.suptitle('PMCA V$_{max}$ rate-limits cytosolic Ca$^{2+}$ recovery '
		f'(EDTA, {scan.length_sec} s)', fontsize=12)
	fig.tight_layout()
	fig.savefig(png_path, dpi=140, bbox_inches='tight')
	plt.close(fig)


def plot_mcu(scan: PerturbationScan, png_path: str) -> None:
	"""2-panel: cyt traces (buffering) + DTS traces (no store rescue)."""
	t = np.arange(scan.cyt.shape[1])
	# Explicit colours: KO=red, baseline=black, over-expression=blue.
	colours = {0.0: '#cc0000', 1.0: '#222222', 4.0: '#1155cc'}
	labels = {0.0: 'MCU knockout (×0)', 1.0: 'baseline (×1)',
		4.0: 'over-expression (×4)'}

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
	for i, f in enumerate(scan.factors):
		c = colours.get(f, None)
		lbl = labels.get(f, f'×{f:g}')
		ax1.plot(t, scan.cyt[i], color=c, linewidth=1.8, label=lbl)
		ax2.plot(t, scan.dts[i], color=c, linewidth=1.8, label=lbl)
	ax1.set_xlabel('Time (s)')
	ax1.set_ylabel(r'Cytosolic Ca$^{2+}$ (nM)')
	ax1.set_title('MCU buffers cytosolic Ca$^{2+}$')
	ax1.legend(fontsize=9)
	ax1.grid(alpha=0.3)

	ax2.set_xlabel('Time (s)')
	ax2.set_ylabel(r'DTS Ca$^{2+}$ (µM)')
	ax2.set_title('DTS depletes regardless of MCU')
	ax2.legend(fontsize=9)
	ax2.grid(alpha=0.3)

	fig.suptitle('MCU buffers cytosolic Ca$^{2+}$ without rescuing the DTS '
		f'(+Ca$^{{2+}}$, {scan.length_sec} s)', fontsize=12)
	fig.tight_layout()
	fig.savefig(png_path, dpi=140, bbox_inches='tight')
	plt.close(fig)


_PLOTTERS = {'pmca': plot_pmca, 'mcu': plot_mcu}


def write_outputs(scan: PerturbationScan, out_path: str) -> None:
	"""Write the NPZ + trace figure for one experiment."""
	scan.to_npz(os.path.join(out_path, f'{scan.exp_key}.npz'))
	_PLOTTERS[scan.exp_key](scan, os.path.join(out_path, f'{scan.exp_key}_traces.png'))


def write_summary(scans: list[PerturbationScan], out_path: str) -> None:
	payload = {
		'experiments': {
			s.exp_key: {
				'title': s.cfg['title'],
				'knob': f"{s.cfg['dict_name']}['{s.cfg['knob_key']}']",
				'baseline_value': s.baseline_value,
				'ca_ex_mM': s.cfg['ca_ex_mM'],
				'length_sec': s.length_sec,
				'conditions': s.scalars,
			}
			for s in scans
		}
	}
	with open(os.path.join(out_path, 'perturbation_summary.json'), 'w') as f:
		json.dump(payload, f, indent='\t')


# ── CLI ─────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		description='PMCA / MCU single-mechanism perturbation experiments (#53).')
	p.add_argument('sim_outdir', nargs='?', default=None,
		help='Output dir under out/. Default = perturbation_<timestamp>.')
	p.add_argument('--experiment', choices=['pmca', 'mcu', 'both'], default='both',
		help='Which experiment(s) to run. Default = both.')
	p.add_argument('--length', '--length-sec', dest='length_sec', type=int,
		default=None, help='Override sim length (s) for all experiments.')
	p.add_argument('--keep-cell-output', action='store_true',
		help='Keep full simOut per condition (default: prune to CalciumTrace).')
	p.add_argument('--no-log-to-shell', dest='log_to_shell',
		action='store_false', help='Suppress per-condition progress prints.')
	p.set_defaults(log_to_shell=True)
	return p


def main(argv: list[str] | None = None) -> None:
	args = _build_parser().parse_args(argv)
	sim_outdir = args.sim_outdir or 'perturbation_{}'.format(
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)

	keys = ['pmca', 'mcu'] if args.experiment == 'both' else [args.experiment]
	print(f'Perturbation experiments: {", ".join(keys)}')
	print(f'  Output: {sim_path}\n')

	scans = []
	for key in keys:
		print(f'[{key}] {EXPERIMENTS[key]["title"]}')
		scan = run_perturbation(sim_path, key,
			length_override=args.length_sec,
			keep_cell_output=args.keep_cell_output,
			log_to_shell=args.log_to_shell)
		write_outputs(scan, sim_path)
		scans.append(scan)
		print()

	write_summary(scans, sim_path)
	print(f'Outputs in: {sim_path}')


if __name__ == '__main__':
	main()
