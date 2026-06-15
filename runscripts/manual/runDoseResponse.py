"""Agonist dose-response — the in-silico predictions for the DTS-depletion experiment.

The recovery-phase analysis raised an unmeasured question: under sustained agonist,
does the platelet dense tubular system (DTS) Ca2+ store deplete completely (model,
~0 uM) or partially (Dolan estimate, 120-180 uM)? The model makes regime-dependent,
falsifiable predictions that an ER-targeted low-affinity GECI experiment could test
(see reports/design/model-driven-experiment-design-2026-06-14.qmd). Running the
sweep (ADP, the reversible/graded agonist) makes them concrete:

  P1  the store ALWAYS transiently empties at the peak (nadir ~0 uM) regardless of
      dose -- so the depletion *nadir* is NOT the dose discriminator
  P2  the SUSTAINED (recovery-phase) [Ca]_dts is dose-graded: it refills at low
      dose and stays empty at high dose, crossing Dolan's 120-180 uM band at an
      intermediate dose (~0.5 uM ADP) -- this is the measurable discriminator
  P3  the autocrine loops are a COMMITMENT SWITCH: with them on (full v0.61) the
      store stays empty at EVERY dose (committed); removing them
      (autocrine_adp_gain=0, cox1_factor=0) restores the graded, reversible
      dose-response -- so the loops convert graded -> all-or-none
  +   PAC-1 (integrin) is dose-graded; peak cytosolic Ca2+ is store-limited
      (dose-insensitive) -- independent per-cell readouts for the model->expt loop

NOTE thrombin is a poor dose-response agonist here: irreversible PAR1 fully cleaves
at any non-zero dose over a long window, so the response is all-or-none (saturated)
-- itself a model prediction (PAR1 is a committed switch). ADP (reversible P2Y1) is
the graded agonist.

This sweeps ONE agonist over a log dose range, in two conditions (full v0.61 vs
autocrine loops off), runs each long enough to see depletion + recovery, and
writes the predicted curves (figure + npz + json).

Usage:
    PYTHONPATH=$PWD python runscripts/manual/runDoseResponse.py [sim_outdir] \\
        [--agonist thrombin|adp] [--grid 7] [--length 350]

Outputs under out/<sim_outdir>/:
    dose_response_<agonist>.png    4-panel predicted curves
    dose_response_<agonist>.npz    doses + per-condition readout arrays
    dose_response_<agonist>.json   metadata + scalars
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from runscripts.manual.runPlateletSim import resolve_sim_path, run_platelet_sim
from reconstruction.platelet.run_config import RunConfig
import wholecell.utils.filepath as fp
from wholecell.io.tablereader import TableReader

# Agonist → (RunConfig peak field, default log dose range, unit label).
AGONISTS = {
	'thrombin': dict(field='thrombin_peak_nM', lo=0.02, hi=5.0, unit='nM',
		label='thrombin', other=dict(adp_peak_uM=0.0, atp_ex_peak_uM=0.0)),
	'adp': dict(field='adp_peak_uM', lo=0.1, hi=30.0, unit='µM',
		label='ADP', other=dict(thrombin_peak_nM=0.0, atp_ex_peak_uM=0.0)),
}

# (key, autocrine_adp_gain, cox1_factor, label, colour).
CONDITIONS = [
	('full', 1.0, 1.0, 'full v0.61 (autocrine loops on)', '#c0392b'),
	('loops_off', 0.0, 0.0, 'autocrine loops off', '#2980b9'),
]

_DOLAN_DTS_BAND_UM = (120.0, 180.0)   # Dolan 2014 Fig 4A dts_min estimate


def _doses(spec: dict, grid: int) -> np.ndarray:
	return np.geomspace(spec['lo'], spec['hi'], grid)


def _harvest(simout_dir: str) -> dict[str, float]:
	cal = TableReader(os.path.join(simout_dir, 'CalciumTrace'))
	dts = cal.readColumn('ca_dts_uM').flatten()
	cyt = cal.readColumn('ca_cyt_nM').flatten()
	ip3 = cal.readColumn('ip3_nM').flatten()
	itg = TableReader(os.path.join(simout_dir, 'IntegrinTrace'))
	pac1 = itg.readColumn('active_frac').flatten()
	tail = max(1, len(dts) // 4)   # last ~quarter = recovery-phase plateau
	return {
		'dts_rest': float(dts[0]),
		'dts_min': float(dts.min()),
		'dts_end': float(dts[-1]),            # P2: recovery-phase refill
		'cyt_peak': float(cyt.max()),
		'ip3_plateau': float(ip3[-tail:].mean()),
		'pac1_max': float(pac1.max()),
	}


def run_dose_response(out_path: str, agonist: str = 'adp', grid: int = 7,
		length: int = 350, log_to_shell: bool = True) -> dict:
	"""Sweep `agonist` over a log dose range in both conditions; harvest readouts."""
	spec = AGONISTS[agonist]
	doses = _doses(spec, grid)
	# results[cond_key][readout] = array over doses
	results: dict[str, dict[str, list]] = {c[0]: {} for c in CONDITIONS}
	keys = ('dts_rest', 'dts_min', 'dts_end', 'cyt_peak', 'ip3_plateau', 'pac1_max')
	for ckey, adp_gain, cox1, _lbl, _col in CONDITIONS:
		rows = []
		for d in doses:
			cell_dir = os.path.join(out_path, f'{ckey}_{agonist}_{d:.4g}')
			fp.makedirs(cell_dir)
			cfg = RunConfig(ca_ex_mM=1.2, autocrine_adp_gain=adp_gain,
				cox1_factor=cox1, **{spec['field']: float(d)}, **spec['other'])
			paths = run_platelet_sim(cell_dir, length_sec=length, seed=0,
				log_to_shell=False, run_config=cfg)
			h = _harvest(paths['sim_out_dir'])
			rows.append(h)
			if log_to_shell:
				print(f'  {ckey:9s} {spec["label"]}={d:7.3g}{spec["unit"]}  '
					f'dts_min={h["dts_min"]:6.1f}  dts_end={h["dts_end"]:6.1f}  '
					f'PAC-1={h["pac1_max"]:.3f}  cyt_pk={h["cyt_peak"]:5.0f}')
		for k in keys:
			results[ckey][k] = [r[k] for r in rows]
	return {'agonist': agonist, 'unit': spec['unit'], 'label': spec['label'],
		'doses': doses.tolist(), 'length_sec': length, 'results': results}


def plot_dose_response(payload: dict, png_path: str) -> None:
	doses = np.array(payload['doses'])
	res = payload['results']
	unit, label = payload['unit'], payload['label']
	full, off = res['full'], res['loops_off']
	C_FULL, C_OFF = '#c0392b', '#2980b9'

	fig, axes = plt.subplots(2, 2, figsize=(12, 9))
	ax_recov, ax_nadir, ax_pac1, ax_cyt = axes.flat

	# P2/P3 (headline) — sustained, recovery-phase store [Ca]_DTS vs dose.
	ax_recov.axhspan(*_DOLAN_DTS_BAND_UM, color='#bdc3c7', alpha=0.45,
		label='Dolan residual band (120–180 µM)')
	ax_recov.plot(doses, off['dts_end'], 's--', color=C_OFF, lw=2,
		label='autocrine loops off (graded, reversible)')
	ax_recov.plot(doses, full['dts_end'], 'o-', color=C_FULL, lw=2,
		label='full v0.61 (committed → stays empty)')
	ax_recov.set_ylabel(r'sustained [Ca²⁺]$_{DTS}$ (recovery, µM)')
	ax_recov.set_title('Sustained store depletion vs dose — the discriminator',
		fontsize=10)
	ax_recov.legend(frameon=False, fontsize=7.5, loc='upper right')

	# P1 — transient nadir at the peak: flat ~0, dose-insensitive.
	ax_nadir.plot(doses, off['dts_min'], 's--', color=C_OFF, lw=2)
	ax_nadir.plot(doses, full['dts_min'], 'o-', color=C_FULL, lw=2)
	ax_nadir.set_ylabel(r'transient store nadir [Ca²⁺]$_{DTS}$ (µM)')
	ax_nadir.set_title('Store nadir at the peak — always ~0 (not the readout)',
		fontsize=10)
	ax_nadir.set_ylim(bottom=-5)

	# PAC-1 (integrin) — dose-graded per-cell readout.
	ax_pac1.plot(doses, off['pac1_max'], 's--', color=C_OFF, lw=2,
		label='loops off')
	ax_pac1.plot(doses, full['pac1_max'], 'o-', color=C_FULL, lw=2,
		label='full v0.61')
	ax_pac1.set_ylabel(r'$\alpha_{IIb}\beta_3$ active fraction (PAC-1$^+$)')
	ax_pac1.set_title('PAC-1 (integrin) dose-response', fontsize=10)
	ax_pac1.legend(frameon=False, fontsize=7.5, loc='lower right')

	# Peak cytosolic Ca²⁺ — store-limited, dose-insensitive.
	ax_cyt.plot(doses, off['cyt_peak'], 's--', color=C_OFF, lw=2)
	ax_cyt.plot(doses, full['cyt_peak'], 'o-', color=C_FULL, lw=2)
	ax_cyt.set_ylabel(r'peak cytosolic Ca²⁺ (nM)')
	ax_cyt.set_title('Peak cytosolic Ca²⁺ — store-limited (flat)', fontsize=10)

	for ax in axes.flat:
		ax.set_xscale('log')
		ax.set_xlabel(f'{label} dose ({unit})')
		ax.grid(alpha=0.3, which='both')

	fig.suptitle(
		f'{label} dose-response — in-silico predictions for the DTS-depletion '
		'experiment', fontsize=13)
	fig.text(0.5, -0.01,
		'Predictions an ER-targeted low-affinity GECI experiment would test. '
		'Top-left (headline): the SUSTAINED store level (recovery phase) is the '
		'measurable discriminator — with the autocrine loops off it grades with '
		'dose and crosses Dolan\'s residual band at ~0.5 µM ADP (graded, '
		'reversible), whereas the full v0.61 model stays empty at EVERY dose: the '
		'autocrine loops are a commitment switch (graded → all-or-none). Top-right: '
		'the transient store nadir is ~0 at all doses (the store always empties at '
		'the peak), so the nadir is not the readout — the recovery is. Bottom-left: '
		'PAC-1 (integrin) is dose-graded (loops add drive). Bottom-right: peak '
		'cytosolic Ca²⁺ is store-limited and dose-insensitive. Aggregation is '
		'inter-cellular and out of scope for a single-cell model.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.03, 1, 0.96])
	fig.savefig(png_path, dpi=150, bbox_inches='tight')
	plt.close(fig)


def write_outputs(payload: dict, out_path: str) -> None:
	agonist = payload['agonist']
	arrays = {'doses': np.array(payload['doses'])}
	for ckey, cond in payload['results'].items():
		for k, v in cond.items():
			arrays[f'{ckey}__{k}'] = np.array(v)
	np.savez(os.path.join(out_path, f'dose_response_{agonist}.npz'), **arrays)
	plot_dose_response(payload,
		os.path.join(out_path, f'dose_response_{agonist}.png'))
	meta = {k: payload[k] for k in ('agonist', 'unit', 'label', 'length_sec')}
	meta['dolan_dts_band_uM'] = list(_DOLAN_DTS_BAND_UM)
	meta['results'] = payload['results']
	meta['doses'] = payload['doses']
	with open(os.path.join(out_path, f'dose_response_{agonist}.json'), 'w') as f:
		json.dump(meta, f, indent='\t')


def _build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(
		description='Agonist dose-response: in-silico predictions for the '
			'DTS-depletion experiment (depth, refill, PAC-1, peak Ca²⁺).')
	p.add_argument('sim_outdir', nargs='?', default=None,
		help='Output dir under out/. Default = dose_response_<timestamp>.')
	p.add_argument('--agonist', choices=sorted(AGONISTS), default='adp',
		help='Agonist to sweep (default adp — the graded, reversible agonist; '
			'thrombin is all-or-none via irreversible PAR1, see module docstring).')
	p.add_argument('--grid', type=int, default=7,
		help='Number of log-spaced doses (default 7).')
	p.add_argument('--length', '--length-sec', dest='length_sec', type=int,
		default=350, help='Simulation length (s). Default 350.')
	p.add_argument('--no-log-to-shell', dest='log_to_shell',
		action='store_false')
	p.set_defaults(log_to_shell=True)
	return p


def main(argv: list[str] | None = None) -> None:
	args = _build_parser().parse_args(argv)
	sim_outdir = args.sim_outdir or 'dose_response_{}'.format(
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)
	print(f'Dose-response: {args.agonist} ({args.grid} doses, {args.length_sec} s)')
	print(f'  Output: {sim_path}\n')
	payload = run_dose_response(sim_path, agonist=args.agonist, grid=args.grid,
		length=args.length_sec, log_to_shell=args.log_to_shell)
	write_outputs(payload, sim_path)
	print(f'\n  Wrote dose_response_{args.agonist}.png / .npz / .json to {sim_path}')


if __name__ == '__main__':
	main()
