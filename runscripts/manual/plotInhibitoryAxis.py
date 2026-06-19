"""Figures for the v0.7 inhibitory axis — P2Y12 / Gi / cAMP / PKA (issue #10).

Two self-contained figures, each running its own simulation(s) with supported
RunConfig knobs (no monkeypatching) and writing a PNG to the output directory
(default ``reports/figures/v0.7``):

  1. mechanism   -> inhibitory_axis_mechanism.png
       Under a standard agonist transient: secreted ADP activates P2Y12, Gi
       lowers cAMP, PKA falls, and the clinical VASP/PRI readout (phospho-VASP)
       falls with it. The new "off" machinery in action.

  2. treatments  -> antiplatelet_treatments.png
       Control vs aspirin (COX-1 = 0) vs clopidogrel (P2Y12 block) vs both,
       under the same activation. PAC-1 integrin activation (the requested
       readout) + cAMP + VASP/PRI + pericellular TXA2, so each drug shows its
       characteristic, biologically-correct signature:
         - clopidogrel: keeps cAMP/VASP-P high -> PKA brake intact -> lower PAC-1
         - aspirin:     abolishes TXA2 (its arm); little PAC-1 effect here

Figure conventions: matplotlib mathtext for chemical formulae (raw unicode
super/subscripts drop glyphs), detailed per-series legends, and a caption
stating the takeaway (standalone thesis artefacts).

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotInhibitoryAxis.py \\
        [--outdir reports/figures/v0.7] [--figure mechanism|treatments|all]
"""

from __future__ import annotations

import argparse
import os
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	_UM_PER_COUNT_EX,
)
from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.io.tablereader import TableReader


def _run(length=250, **config_kwargs):
	"""Run one platelet sim into a temp dir; return its simOut path."""
	run_config = RunConfig(ca_ex_mM=1.2, **config_kwargs)
	td = tempfile.mkdtemp()
	paths = run_platelet_sim(td, length_sec=length, seed=0,
		log_to_shell=False, run_config=run_config)
	return paths['sim_out_dir']


def _trace(sim_out, column):
	return TableReader(os.path.join(sim_out, 'CalciumTrace')).readColumn(
		column).flatten()


def _bulk(sim_out):
	rb = TableReader(os.path.join(sim_out, 'BulkMolecules'))
	return list(rb.readAttribute('objectNames')), rb.readColumn('counts')


def _pac1(sim_out):
	"""αIIbβ3 active fraction (%) — the per-cell PAC-1 readout."""
	ids, counts = _bulk(sim_out)
	act = counts[:, ids.index('aIIbb3_active[pl]')].astype(float)
	rest = counts[:, ids.index('aIIbb3_resting[pl]')].astype(float)
	total = act + rest
	return 100.0 * np.divide(act, total, out=np.zeros_like(act), where=total > 0)


def _txa2_uM(sim_out):
	"""Pericellular synthesised TXA2 (µM) — the aspirin (COX-1) arm."""
	ids, counts = _bulk(sim_out)
	return counts[:, ids.index('TXA2[e]')].astype(float) * _UM_PER_COUNT_EX


# ── 1. Inhibitory-axis mechanism ─────────────────────────────────────────────

def fig_mechanism(outdir: str) -> str:
	"""The cAMP/PKA off-switch responding to agonist (standard transient, WT)."""
	sim = _run(250)
	t = np.arange(len(_trace(sim, 'camp_uM')))
	camp = _trace(sim, 'camp_uM')
	pka = _trace(sim, 'pka_active_frac')
	vaspp = _trace(sim, 'vasp_phos_frac')
	p2y12 = _trace(sim, 'p2y12_active_frac')

	fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.7))

	# Left: P2Y12 occupancy drives cAMP down.
	a1.plot(t, p2y12, lw=2.2, color='#c0392b',
		label='P2Y12 active fraction (ADP-bound)')
	a1.set_xlabel('time (s)')
	a1.set_ylabel('P2Y12 active fraction', color='#c0392b')
	a1.tick_params(axis='y', labelcolor='#c0392b')
	a1.set_ylim(0, 1)
	a1r = a1.twinx()
	a1r.plot(t, camp, lw=2.4, color='#2980b9', label=r'cAMP ($\mu$M)')
	a1r.axhline(1.0, color='#2980b9', lw=0.8, ls=':', alpha=0.6)
	a1r.set_ylabel(r'cAMP ($\mu$M)', color='#2980b9')
	a1r.tick_params(axis='y', labelcolor='#2980b9')
	a1r.set_ylim(bottom=0)
	a1.set_title(r'ADP $\rightarrow$ P2Y12 $\rightarrow$ G$_i$ lowers cAMP',
		fontsize=10)
	la, lba = a1.get_legend_handles_labels()
	lb, lbb = a1r.get_legend_handles_labels()
	a1.legend(la + lb, lba + lbb, loc='center right', fontsize=8, frameon=False)
	a1.grid(alpha=0.3)

	# Right: PKA and the VASP/PRI readout fall together.
	a2.plot(t, pka, lw=2.4, color='#8e44ad', label='PKA active fraction')
	a2.plot(t, vaspp, lw=2.4, color='#27ae60',
		label='phospho-VASP fraction (VASP/PRI)')
	a2.axhline(vaspp[0], color='#27ae60', lw=0.8, ls=':', alpha=0.6,
		label=f'resting VASP-P = {vaspp[0]:.2f}')
	a2.set_xlabel('time (s)')
	a2.set_ylabel('active / phosphorylated fraction')
	a2.set_ylim(0, 1)
	a2.set_title(r'Falling cAMP $\rightarrow$ less PKA $\rightarrow$ less VASP-P',
		fontsize=10)
	a2.legend(loc='center right', fontsize=8, frameon=False)
	a2.grid(alpha=0.3)

	fig.suptitle('v0.7 inhibitory axis — P2Y12 / G$_i$ / cAMP / PKA '
		r'(standard thrombin + ADP + ATP transient, +$\mathrm{Ca^{2+}}$ 1.2 mM)',
		fontsize=12)
	fig.text(0.5, -0.02,
		'Secreted and applied ADP binds P2Y12 (left, red); G$_i$ inhibits '
		'adenylyl cyclase so cAMP falls from its ~1 µM resting tone (left, blue). '
		'Lower cAMP means less PKA activity and, with it, less phospho-VASP '
		'(right) — the basis of the clinical VASP/PRI assay for P2Y12 inhibition. '
		'At rest (ADP = 0) P2Y12 is silent, cAMP sits at basal, and the PKA brake '
		'is exactly neutral, so the resting fixed point and Dolan 5/5 are '
		'unchanged.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.02, 1, 0.95])
	out = os.path.join(outdir, 'inhibitory_axis_mechanism.png')
	fig.savefig(out, dpi=150, bbox_inches='tight')
	plt.close(fig)
	return out


# ── 2. Antiplatelet treatments ───────────────────────────────────────────────

_TREATMENTS: list[tuple[str, dict[str, float], str, str]] = [
	('control',             {},                                       '#2c3e50', '-'),
	('aspirin (COX-1 = 0)', dict(cox1_factor=0.0),                    '#e67e22', '--'),
	('clopidogrel (P2Y12)', dict(p2y12_block=1.0),                    '#c0392b', '-'),
	('aspirin + clopidogrel', dict(cox1_factor=0.0, p2y12_block=1.0), '#8e44ad', ':'),
]


def fig_treatments(outdir: str) -> str:
	"""Control vs aspirin vs clopidogrel vs both — PAC-1 + cAMP + VASP + TXA2."""
	runs = [(label, _run(250, **kw), color, ls)
		for label, kw, color, ls in _TREATMENTS]
	t = np.arange(len(_pac1(runs[0][1])))

	fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))
	(axA, axB), (axC, axD) = axes

	for label, sim, color, ls in runs:
		axA.plot(t, _pac1(sim), lw=2.2, color=color, ls=ls, label=label)
		axB.plot(t, _trace(sim, 'camp_uM'), lw=2.2, color=color, ls=ls, label=label)
		axC.plot(t, _trace(sim, 'vasp_phos_frac'), lw=2.2, color=color, ls=ls, label=label)
		axD.plot(t, _txa2_uM(sim), lw=2.2, color=color, ls=ls, label=label)

	axA.set_title('A — PAC-1 integrin (αIIbβ3) activation', fontsize=10, fontweight='bold')
	axA.set_ylabel('active αIIbβ3 (% — PAC-1 readout)')
	axA.set_ylim(bottom=0)

	axB.set_title('B — cAMP', fontsize=10, fontweight='bold')
	axB.set_ylabel(r'cAMP ($\mu$M)')
	axB.set_ylim(bottom=0)

	axC.set_title('C — phospho-VASP (VASP/PRI assay)', fontsize=10, fontweight='bold')
	axC.set_ylabel('phospho-VASP fraction')
	axC.set_ylim(0, 1)

	axD.set_title(r'D — pericellular $\mathrm{TXA_2}$ (aspirin target)', fontsize=10, fontweight='bold')
	axD.set_ylabel(r'$[\mathrm{TXA_2}]_\mathrm{e}$ ($\mu$M)')
	axD.set_ylim(bottom=0)

	for ax in (axA, axB, axC, axD):
		ax.set_xlabel('time (s)')
		ax.grid(alpha=0.3)
		ax.legend(loc='best', fontsize=8, frameon=False)

	fig.suptitle('Antiplatelet treatments in the platelet whole-cell model '
		r'(standard agonist, +$\mathrm{Ca^{2+}}$ 1.2 mM)',
		fontsize=13, fontweight='bold')
	fig.text(0.5, -0.01,
		'Each drug shows its mechanism-specific signature. Clopidogrel '
		'(P2Y12 antagonist) blocks ADP-driven cAMP lowering, so cAMP (B) and '
		'phospho-VASP (C) stay at their high resting tone — the PKA brake remains '
		'engaged and PAC-1 integrin activation (A) is reduced versus control: the '
		'real clopidogrel mechanism and the basis of the VASP/PRI assay. Aspirin '
		'(COX-1 knockout) abolishes thromboxane synthesis (D, $\\mathrm{TXA_2}$ '
		'$\\rightarrow$ 0) but has little effect on PAC-1 here, because the '
		'$\\mathrm{TXA_2}\\rightarrow$TP$\\rightarrow$G$_q$ amplification it '
		'removes is secondary to the primary agonist drive in this regime. The '
		'cytosolic-$\\mathrm{Ca^{2+}}$ effect of P2Y12 is small (the '
		'store-limited/SOCE clamp), '
		'so the visible P2Y12 action is on the downstream integrin output.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.03, 1, 0.96])
	out = os.path.join(outdir, 'antiplatelet_treatments.png')
	fig.savefig(out, dpi=150, bbox_inches='tight')
	plt.close(fig)
	return out


def main(argv: list[str] | None = None) -> None:
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--outdir', default='reports/figures/v0.7')
	parser.add_argument('--figure', default='all',
		choices=['mechanism', 'treatments', 'all'])
	args = parser.parse_args(argv)

	os.makedirs(args.outdir, exist_ok=True)
	if args.figure in ('mechanism', 'all'):
		print('Wrote', fig_mechanism(args.outdir))
	if args.figure in ('treatments', 'all'):
		print('Wrote', fig_treatments(args.outdir))


if __name__ == '__main__':
	main()
