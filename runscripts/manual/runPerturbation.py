"""Single-mechanism perturbation experiments (PMCA V_max, MCU uptake, PKC).

In-silico knockdown/knockout scans on the calcium model, producing the
mechanistic-finding figures the thesis Results section claims (issue #53;
the PKC scan is the v0.6 #57 deliverable). Each condition overrides one
rate constant in `calcium_signalling`, runs a platelet sim, harvests the
cytosolic / DTS Ca²⁺ traces (plus an optional experiment-specific aux
column), then restores the constant — the same live-override pattern
`runDoseSweep.py` uses.

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
        (full in every case — MCU doesn't refill the store). NB MCU=0 *raises*
        cyt in this buffer-only model, whereas real MCU⁻/⁻ platelets show
        *reduced* agonist-evoked cyt Ca²⁺ (Ghatge 2026; Ajanel 2025) — the
        model diverges; see issue #76.

Exp C — PKC desensitises P2Y1, throttling the ADP arm (v0.6, issue #57)
    Knob   : K_P2Y1_DES['k_des']  (P2Y1 desensitisation rate; 0 = knockout)
    Scan   : {0 (KO), 1 (baseline feedback)}
    Protocol: +Ca²⁺, ADP-only (thrombin off) to isolate the P2Y1 arm — the
        thrombin-driven PARs otherwise dominate Gαq and mask the effect.
    Observable: P2Y1 desensitised fraction (aux column) — stays ~0 in the
        knockout, rises to ~0.5–0.7 with feedback — plus cytosolic Ca²⁺. The
        knockout reproduces the increased-reactivity phenotype of
        desensitisation-resistant P2Y1 (Nicholas 2023). The single-transient
        Ca²⁺ amplitude effect is modest (the response is store-limited); the
        receptor desensitisation itself is the clear, measured readout.

Exp D — PKC → PLCβ phosphorylation damps IP₃ (v0.6 Slice 3, Purvis route)
    Knob   : K_PLCB_PHOS['k_plcb_phos']  (PLCβ phosphorylation; 0 = knockout)
    Scan   : {0 (KO), 1 (baseline feedback)}
    Protocol: +Ca²⁺, standard thrombin + ADP. Unlike P2Y1 desensitisation,
        this brake sits on the SHARED PLCβ node, so it damps the cascade
        downstream of all receptors and the standard transient shows it.
    Observable: IP₃ (aux) — the knockout plateaus high; the baseline feedback
        peaks then returns toward baseline (Purvis 2008 Fig 2C/5E) — plus the
        PLCβ phosphorylated fraction (aux), which rises to ~0.3 with feedback.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/runPerturbation.py [sim_outdir] \\
        [--experiment pmca|mcu|pkc|both] [--length N] [--keep-cell-output]

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
from typing import Any

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from runscripts.manual.runPlateletSim import resolve_sim_path, run_platelet_sim
from reconstruction.platelet.run_config import RunConfig
from reconstruction.platelet.dataclasses.process import calcium_signalling as cs_mod
import wholecell.utils.filepath as fp
from wholecell.io.tablereader import TableReader


# ── Experiment registry ─────────────────────────────────────────────────────
# Each experiment scales ONE constant via a RunConfig field (config_field): the
# baseline of cs_mod.<dict_name>[<knob_key>] is multiplied by each factor in turn
# inside the ODE (no global mutation). ca_ex_mM and length_sec set the protocol.
EXPERIMENTS: dict[str, dict] = {
	'pmca': dict(
		title='PMCA V_max rate-limits cytosolic Ca recovery (EDTA transient)',
		dict_name='K_PMCA', knob_key='k_cat', config_field='pmca_kcat_scale',
		knob_label='PMCA V_max',
		factors=[0.25, 0.5, 1.0, 2.0, 4.0],
		ca_ex_mM=0.0, length_sec=300,
	),
	'mcu': dict(
		title='MCU buffers cytosolic Ca without rescuing the DTS store (+Ca)',
		dict_name='K_MITO', knob_key='V_max_MCU', config_field='mcu_vmax_scale',
		knob_label='MCU V_max',
		factors=[0.0, 1.0, 4.0],
		ca_ex_mM=1.2, length_sec=400,
	),
	'pkc': dict(
		# v0.6 — PKC-mediated P2Y1 desensitisation (DAG → PKC → P2Y1).
		# Knock out the desensitisation rate (k_des → 0) vs the baseline
		# feedback. ADP-only (thrombin off) isolates the P2Y1 arm, since
		# the thrombin-driven PARs otherwise dominate Gαq and mask the
		# P2Y1-specific effect. Prediction (Nicholas 2023): the knockout
		# keeps P2Y1 active (desensitised fraction stays ~0) and gives a
		# larger / more sustained response than the baseline feedback.
		title='PKC desensitisation of P2Y1 throttles the ADP response '
			'(ADP-only, +Ca)',
		dict_name='K_P2Y1_DES', knob_key='k_des', config_field='k_des_scale',
		knob_label='P2Y1 k_des',
		factors=[0.0, 1.0],
		ca_ex_mM=1.2, length_sec=300,
		sim_kwargs=dict(thrombin_peak_nM=0.0),     # ADP-only: isolate P2Y1
		aux_cols=['p2y1_desensitised_frac'],       # mechanism readout
	),
	'plcb': dict(
		# v0.6 Slice 3 — PKC → PLCβ phosphorylation (Purvis 2008 route).
		# Knock out the phosphorylation rate (k_plcb_phos → 0) vs baseline.
		# Unlike the P2Y1 route this brake sits on the SHARED PLCβ node, so
		# the standard thrombin + ADP transient shows it: PKC sequesters
		# PLCβ out of the Gq-activatable pool, lowering IP₃ — the Purvis
		# "return toward baseline" result. IP₃ is the clear readout (the
		# cytosolic Ca²⁺ amplitude stays store-limited).
		title='PKC → PLCβ phosphorylation damps IP₃ (Purvis route, +Ca)',
		dict_name='K_PLCB_PHOS', knob_key='k_plcb_phos',
		config_field='k_plcb_phos_scale', knob_label='PLCβ k_phos',
		factors=[0.0, 1.0],
		ca_ex_mM=1.2, length_sec=300,
		aux_cols=['ip3_nM', 'plcb_phosphorylated_frac'],
	),
}


def harvest_traces(simout_dir: str) -> tuple[np.ndarray, np.ndarray]:
	"""Read the cytosolic (nM) and DTS (µM) Ca²⁺ time series from CalciumTrace."""
	reader = TableReader(os.path.join(simout_dir, 'CalciumTrace'))
	cyt = reader.readColumn('ca_cyt_nM').flatten()
	dts = reader.readColumn('ca_dts_uM').flatten()
	return cyt, dts


def harvest_column(simout_dir: str, col: str) -> np.ndarray:
	"""Read one named CalciumTrace column (e.g. an experiment's aux readout)."""
	reader = TableReader(os.path.join(simout_dir, 'CalciumTrace'))
	return reader.readColumn(col).flatten()


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
	# Optional experiment-specific traces, keyed by CalciumTrace column name;
	# each value is (n_factor, n_time). Saved into the NPZ as `aux_<col>`.
	aux: dict[str, np.ndarray] | None = None

	def to_npz(self, path: str) -> None:
		arrays: dict[str, Any] = dict(
			factors=np.array(self.factors), cyt=self.cyt, dts=self.dts,
			baseline_value=self.baseline_value, length_sec=self.length_sec)
		for col, mat in (self.aux or {}).items():
			arrays[f'aux_{col}'] = mat
		np.savez(path, **arrays)


# ── Driver ────────────────────────────────────────────────────────────────

def run_perturbation(out_path: str, exp_key: str,
		length_override: int | None = None,
		keep_cell_output: bool = False,
		log_to_shell: bool = True) -> PerturbationScan:
	"""Run one experiment's factor scan; always restore the knob afterwards."""
	cfg = EXPERIMENTS[exp_key]
	# Baseline value of the perturbed constant — read-only, for axis labels;
	# the perturbation itself is applied through the RunConfig scale field (the
	# factor multiplies the baseline inside the ODE), not by mutating the dict.
	baseline = float(getattr(cs_mod, cfg['dict_name'])[cfg['knob_key']])
	config_field = cfg['config_field']
	length = length_override or cfg['length_sec']
	factors = list(cfg['factors'])

	cells_root = os.path.join(out_path, f'{exp_key}_cells')
	fp.makedirs(cells_root)

	extra_kwargs = cfg.get('sim_kwargs', {})       # per-experiment agonist overrides
	aux_cols = cfg.get('aux_cols', [])             # optional extra trace columns

	cyt_rows: list[np.ndarray] = []
	dts_rows: list[np.ndarray] = []
	aux_rows: dict[str, list[np.ndarray]] = {c: [] for c in aux_cols}
	scalars: list[dict] = []
	for f in factors:
		run_config = RunConfig(
			ca_ex_mM=cfg['ca_ex_mM'], **extra_kwargs, **{config_field: float(f)})
		cell_dir = os.path.join(cells_root, f'x{f:g}')
		fp.makedirs(cell_dir)
		paths = run_platelet_sim(cell_dir, length_sec=length, seed=0,
			log_to_shell=False, run_config=run_config)
		cyt, dts = harvest_traces(paths['sim_out_dir'])
		cyt_rows.append(cyt)
		dts_rows.append(dts)
		sc = {'factor': float(f), **reduce_scalars(cyt, dts)}
		for col in aux_cols:
			trace = harvest_column(paths['sim_out_dir'], col)
			aux_rows[col].append(trace)
			sc[f'{col}_max'] = float(trace.max())
		scalars.append(sc)
		if log_to_shell:
			print(f'  {cfg["knob_label"]} ×{f:<5g} → '
				f'peak cyt {sc["peak_cyt_nM"]:>6.1f} nM, '
				f'recovery-AUC {sc["recovery_tail_auc_nMs"]:>9.0f} nM·s, '
				f'DTS min {sc["dts_min_uM"]:>6.1f} µM')
		if not keep_cell_output:
			_prune_cell(cell_dir)

	return PerturbationScan(
		exp_key=exp_key, cfg=cfg, baseline_value=baseline, length_sec=length,
		factors=factors, cyt=np.array(cyt_rows), dts=np.array(dts_rows),
		scalars=scalars,
		aux={c: np.array(rows) for c, rows in aux_rows.items()} or None)


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


def plot_pkc(scan: PerturbationScan, png_path: str) -> None:
	"""2-panel: cytosolic Ca²⁺ + P2Y1 desensitised fraction, knockout vs baseline."""
	t = np.arange(scan.cyt.shape[1])
	# k_des ×0 = PKC-desensitisation knockout; ×1 = baseline feedback.
	colours = {0.0: '#cc0000', 1.0: '#222222'}
	labels = {0.0: 'PKC knockout (k$_{des}$×0)', 1.0: 'baseline feedback (×1)'}

	des = (scan.aux or {}).get('p2y1_desensitised_frac')
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
	for i, f in enumerate(scan.factors):
		c = colours.get(f, None)
		lbl = labels.get(f, f'×{f:g}')
		ax1.plot(t, scan.cyt[i], color=c, linewidth=1.8, label=lbl)
		if des is not None:
			ax2.plot(t, des[i], color=c, linewidth=1.8, label=lbl)
	ax1.set_xlabel('Time (s)')
	ax1.set_ylabel(r'Cytosolic Ca$^{2+}$ (nM)')
	ax1.set_title('ADP-evoked cytosolic Ca$^{2+}$ (store-limited)')
	ax1.legend(fontsize=9)
	ax1.grid(alpha=0.3)

	ax2.set_xlabel('Time (s)')
	ax2.set_ylabel('P2Y1 desensitised fraction')
	ax2.set_title('PKC desensitises the active P2Y1 receptor')
	ax2.set_ylim(-0.02, 1.0)
	ax2.legend(fontsize=9)
	ax2.grid(alpha=0.3)

	fig.suptitle('PKC-mediated P2Y1 desensitisation: knockout vs baseline '
		f'(ADP-only, +Ca$^{{2+}}$, {scan.length_sec} s)', fontsize=12)
	fig.tight_layout()
	fig.savefig(png_path, dpi=140, bbox_inches='tight')
	plt.close(fig)


def plot_plcb(scan: PerturbationScan, png_path: str) -> None:
	"""2-panel: IP₃ + PLCβ phosphorylated fraction, knockout vs baseline."""
	t = np.arange(scan.cyt.shape[1])
	colours = {0.0: '#cc0000', 1.0: '#222222'}
	labels = {0.0: 'PLCβ-phos knockout (k$_{phos}$×0)',
		1.0: 'baseline feedback (×1)'}
	ip3 = (scan.aux or {}).get('ip3_nM')
	phos = (scan.aux or {}).get('plcb_phosphorylated_frac')

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
	for i, f in enumerate(scan.factors):
		c = colours.get(f, None)
		lbl = labels.get(f, f'×{f:g}')
		if ip3 is not None:
			ax1.plot(t, ip3[i], color=c, linewidth=1.8, label=lbl)
		if phos is not None:
			ax2.plot(t, phos[i], color=c, linewidth=1.8, label=lbl)
	ax1.set_xlabel('Time (s)')
	ax1.set_ylabel('IP$_3$ (nM)')
	ax1.set_title('PKC → PLCβ phosphorylation lowers IP$_3$ (Purvis route)')
	ax1.legend(fontsize=9)
	ax1.grid(alpha=0.3)

	ax2.set_xlabel('Time (s)')
	ax2.set_ylabel('PLCβ phosphorylated fraction')
	ax2.set_title('PKC sequesters PLCβ out of the Gq-activatable pool')
	ax2.set_ylim(-0.02, 1.0)
	ax2.legend(fontsize=9)
	ax2.grid(alpha=0.3)

	fig.suptitle('PKC → PLCβ phosphorylation: knockout vs baseline '
		f'(thrombin + ADP, +Ca$^{{2+}}$, {scan.length_sec} s)', fontsize=12)
	fig.tight_layout()
	fig.savefig(png_path, dpi=140, bbox_inches='tight')
	plt.close(fig)


_PLOTTERS = {'pmca': plot_pmca, 'mcu': plot_mcu, 'pkc': plot_pkc,
	'plcb': plot_plcb}


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
	p.add_argument('--experiment',
		choices=['pmca', 'mcu', 'pkc', 'plcb', 'both'], default='both',
		help='Which experiment(s) to run. Default = both (pmca + mcu); '
			'pkc / plcb are the v0.6 PKC-feedback knockouts (P2Y1 '
			'desensitisation / PLCβ phosphorylation).')
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
