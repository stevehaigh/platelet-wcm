"""DTS-realism figures (v0.7): the V_IM thermodynamic-floor fix + the γ_IP3R
resting-balance recalibration, and what they do (or don't) change.

This version made the dense-tubular-system (DTS) Ca²⁺ store behave honestly:

  1. V_IM = 0 mV  — the ER/DTS membrane holds ≈0 mV (it's counter-ion
     permeable), so the passive IP3R can only equilibrate the store toward
     *cytosolic* free Ca²⁺, never below it. The old V_IM = −60 mV (a
     plasma-membrane value) drove the store sub-cytosolically — impossible
     for a passive channel.
  2. γ_IP3R 0.075 → 0.135 pS — re-derived for V_IM = 0 so an unstimulated
     cell holds a STABLE resting fixed point (~250 µM / ~100 nM) instead of
     drifting. (The V_IM change alone unbalanced the resting IP3R-leak ⇌ SERCA
     coupling, over-filling the store.)

Conclusion (vindicates "View B"): once homeostasis is correct, the store
deeply depletes under sustained saturating agonist — but stays above cytosol —
and a held residual is unachievable via IP3R inactivation. Partial depletion
appears at lower / threshold doses (the dose-response figure).

V_IM and γ_IP3R are module constants (not RunConfig fields), so this figure
tool sets them around each run with save/restore — the only place we vary them.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotDtsRealism.py
        [--out-dir reports/figures/v0.7]
"""

from __future__ import annotations

import argparse
import os
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import reconstruction.platelet.dataclasses.process.calcium_signalling as cs
import wholecell.utils.filepath as fp
from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.io.tablereader import TableReader

# Parameter sets (V_IM in mV, γ in pS). NEW = the committed v0.7 defaults.
OLD       = (-60.0, 0.075)   # pre-session main (V_IM=−60, γ=0.075)
VIM0_ONLY = (0.0,   0.075)   # V_IM=0 but γ not yet re-derived (broken rest)
NEW       = (0.0,   0.135)   # v0.7: V_IM=0 + recalibrated γ

_COLS = ('ca_cyt_nM', 'ca_dts_uM', 'ip3_nM', 'soce_flux_nMs')


def run_with(v_im_mV: float, gamma_pS: float, ca_ex_mM: float = 1.2,
		length: int = 200, at_rest: bool = False,
		adp_peak_uM: float | None = None,
		thrombin_off: bool = False,
		loops_off: bool = False) -> dict[str, np.ndarray]:
	"""Run one sim with V_IM/γ patched in (save/restore); harvest traces."""
	old_vim, old_g = cs.V_IM_V, cs.GAMMA_IP3R_S
	cs.V_IM_V = v_im_mV / 1000.0
	cs.GAMMA_IP3R_S = gamma_pS * 1e-12
	try:
		kw: dict = dict(ca_ex_mM=ca_ex_mM)
		if at_rest:
			kw.update(thrombin_peak_nM=0.0, adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
		if adp_peak_uM is not None:
			kw['adp_peak_uM'] = adp_peak_uM
		if thrombin_off:
			kw.update(thrombin_peak_nM=0.0, atp_ex_peak_uM=0.0)
		if loops_off:
			kw.update(autocrine_adp_gain=0.0, cox1_factor=0.0)
		td = tempfile.mkdtemp()
		paths = run_platelet_sim(td, length_sec=length, seed=0,
			log_to_shell=False, run_config=RunConfig(**kw))
		r = TableReader(os.path.join(paths['sim_out_dir'], 'CalciumTrace'))
		return {c: r.readColumn(c).flatten() for c in _COLS}
	finally:
		cs.V_IM_V, cs.GAMMA_IP3R_S = old_vim, old_g


# ── Figure 1: the thermodynamic floor (V_IM) ──────────────────────────────

def fig_floor(out_dir: str) -> None:
	"""At the recalibrated γ, vary V_IM: −60 mV drives the store below
	cytosol; 0 mV bottoms it at the cytosolic equilibrium."""
	new = run_with(0.0, NEW[1], length=200)            # V_IM = 0
	oldv = run_with(-60.0, NEW[1], length=200)         # V_IM = −60 (same γ)
	t = np.arange(len(new['ca_dts_uM']))
	fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

	ax = axes[0]
	ax.semilogy(t, oldv['ca_dts_uM'], '#c0392b', lw=2,
		label=r'$V_{IM}=-60$ mV (old) — drains below cytosol')
	ax.semilogy(t, new['ca_dts_uM'], '#2980b9', lw=2,
		label=r'$V_{IM}=0$ mV (v0.7) — bottoms at equilibrium')
	cyt_floor = new['ca_cyt_nM'].max() / 1000.0
	ax.axhline(cyt_floor, color='#7f8c8d', ls=':', lw=1.3,
		label=r'cytosolic free $\mathrm{Ca^{2+}}$ (the thermodynamic floor)')
	ax.set_xlabel('time (s)')
	ax.set_ylabel(r'DTS store free $\mathrm{Ca^{2+}}$ ($\mu$M, log)')
	ax.set_title(r'DTS store depletion', fontsize=10)
	ax.legend(frameon=False, fontsize=7.5, loc='upper right')
	ax.grid(alpha=0.3, which='both')

	ax = axes[1]
	ax.plot(t, oldv['ca_cyt_nM'], '#c0392b', lw=2, label=r'$V_{IM}=-60$ mV')
	ax.plot(t, new['ca_cyt_nM'], '#2980b9', lw=2, label=r'$V_{IM}=0$ mV')
	ax.set_xlabel('time (s)')
	ax.set_ylabel(r'cytosolic $\mathrm{Ca^{2+}}$ (nM)')
	ax.set_title(r'Cytosolic $\mathrm{Ca^{2+}}$ — essentially unchanged', fontsize=10)
	ax.legend(frameon=False, fontsize=8, loc='upper right')
	ax.grid(alpha=0.3)

	fig.suptitle(r'The IP$_3$R is a passive channel: with $V_{IM}=0$ the DTS '
		r'store can no longer drain below cytosolic free $\mathrm{Ca^{2+}}$',
		fontsize=12)
	fig.text(0.5, -0.05,
		'The IP3R Ca2+ flux is gradient-driven (Nernst), so at zero net flux the '
		'store sits at [Ca]_DTS/[Ca]_cyt = exp(V_IM·zF/RT). The old V_IM = −60 mV '
		'(a plasma-membrane value borrowed from Dolan) made that ratio ~0.01, '
		'actively pulling the free store to ~0.03 µM — BELOW the ~0.4 µM cytosol '
		'(left, red), which a passive channel cannot physically do. Setting '
		'V_IM = 0 mV (the ER/DTS membrane is counter-ion permeable and holds ≈0 mV) '
		'makes the equilibrium [Ca]_DTS = [Ca]_cyt, so the store bottoms at the '
		'cytosolic floor (~1 µM, SERCA holding it just above; left, blue). The '
		'cytosolic transient is essentially unchanged (right) — peak is clamped by '
		'buffering / PMCA / SOCE, so Dolan 5/5 is unaffected. Same recalibrated '
		'γ_IP3R = 0.135 pS in both traces; only V_IM differs.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.05, 1, 0.95])
	path = os.path.join(out_dir, 'dts_thermodynamic_floor.png')
	fig.savefig(path, dpi=150, bbox_inches='tight')
	plt.close(fig)
	print(f'  wrote {path}')


# ── Figure 2: resting homeostasis (γ recalibration) ───────────────────────

def fig_resting(out_dir: str) -> None:
	"""Unstimulated 400 s: old drifts, V_IM=0-only drifts worse, recalibrated
	γ holds a stable fixed point."""
	L = 400
	old = run_with(*OLD, length=L, at_rest=True)
	vim0 = run_with(*VIM0_ONLY, length=L, at_rest=True)
	new = run_with(*NEW, length=L, at_rest=True)
	t = np.arange(len(new['ca_dts_uM']))
	series = [
		(old,  '#c0392b', '-',
			r'old: $V_{IM}=-60$, $\gamma=0.075$ pS'),
		(vim0, '#e67e22', '--',
			r'$V_{IM}=0$ only, $\gamma=0.075$ pS (rest unbalanced)'),
		(new,  '#2980b9', '-',
			r'v0.7: $V_{IM}=0$, $\gamma=0.135$ pS (recalibrated)'),
	]
	fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
	for d, c, ls, lab in series:
		axes[0].plot(t, d['ca_dts_uM'], color=c, ls=ls, lw=2, label=lab)
		axes[1].plot(t, d['ca_cyt_nM'], color=c, ls=ls, lw=2, label=lab)
	axes[0].axhspan(100, 400, color='#bdc3c7', alpha=0.30,
		label=r'Dolan resting band (100–400 $\mu$M)')
	axes[0].set_ylabel(r'DTS store $\mathrm{Ca^{2+}}$ ($\mu$M)')
	axes[0].set_title('DTS store at rest', fontsize=10)
	axes[1].axhspan(25, 100, color='#bdc3c7', alpha=0.30,
		label='Dolan resting band (25–100 nM)')
	axes[1].set_ylabel(r'cytosolic $\mathrm{Ca^{2+}}$ (nM)')
	axes[1].set_title(r'Cytosolic $\mathrm{Ca^{2+}}$ at rest', fontsize=10)
	for ax in axes:
		ax.set_xlabel('time (s)')
		ax.legend(frameon=False, fontsize=7.5, loc='best')
		ax.grid(alpha=0.3)
	fig.suptitle('Resting homeostasis: re-deriving γ_IP3R for V_IM=0 restores a '
		'stable resting fixed point', fontsize=12)
	fig.text(0.5, -0.05,
		'An unstimulated cell should hold ~250 µM DTS / ~100 nM cytosol. The store '
		'balance is the resting IP3R leak ⇌ SERCA; γ_IP3R was originally tuned for '
		'V_IM = −60 mV. Moving to V_IM = 0 mV shrinks the resting driving force, so '
		'the old γ = 0.075 pS under-leaks and SERCA over-fills the store — the cell '
		'drifts to ~500 µM DTS / ~6 nM cytosol (orange, worse than the pre-session '
		'baseline in red). Re-deriving γ to 0.135 pS (still within the 0.05–0.5 pS '
		'plausibility range) restores a stable fixed point: DTS 250→255 µM, cytosol '
		'100→96 nM, settled over 400 s (blue). This pre-existing resting drift was '
		'masked because validation runs initialise at 250/100 and run short, '
		'agonist-driven windows.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.05, 1, 0.95])
	path = os.path.join(out_dir, 'resting_homeostasis.png')
	fig.savefig(path, dpi=150, bbox_inches='tight')
	plt.close(fig)
	print(f'  wrote {path}')


# ── Figure 3: calibration sweeps ──────────────────────────────────────────

def fig_sweeps(out_dir: str) -> None:
	"""V_IM sets the floor; γ trades resting stability against active residual."""
	# Panel A: V_IM sweep (at new γ) → DTS min (+Ca).
	vims = [-60, -40, -20, 0, 40]
	dts_min, cyt_peak = [], []
	for v in vims:
		d = run_with(float(v), NEW[1], length=200)
		dts_min.append(d['ca_dts_uM'].min())
		cyt_peak.append(d['ca_cyt_nM'].max())
	# Panel B: γ sweep (at V_IM=0) → resting DTS (stability) + active residual.
	gammas = [0.075, 0.10, 0.135, 0.17, 0.22]
	rest_dts, resid = [], []
	for g in gammas:
		rest_dts.append(run_with(0.0, g, length=400, at_rest=True)['ca_dts_uM'][-1])
		resid.append(run_with(0.0, g, length=200)['ca_dts_uM'][-1])

	fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
	ax = axes[0]
	ax.semilogy(vims, dts_min, 'o-', color='#2980b9', lw=2, ms=7,
		label=r'DTS min (+Ca)')
	ax.axhline(cyt_peak[0] / 1000.0, color='#7f8c8d', ls=':', lw=1.3,
		label=r'cytosolic floor')
	ax.axvline(0, color='#27ae60', ls='--', lw=1.3, label=r'v0.7 ($V_{IM}=0$)')
	ax.set_xlabel(r'$V_{IM}$ (mV)')
	ax.set_ylabel(r'DTS store minimum ($\mu$M, log)')
	ax.set_title(r'$V_{IM}$ sets the depletion FLOOR (all $V_{IM}$: Dolan 5/5)',
		fontsize=9.5)
	ax.legend(frameon=False, fontsize=8, loc='upper left')
	ax.grid(alpha=0.3, which='both')

	ax = axes[1]
	ax.plot(gammas, rest_dts, 's-', color='#c0392b', lw=2, ms=7,
		label='resting DTS (400 s) — stability')
	ax.plot(gammas, resid, 'o-', color='#8e44ad', lw=2, ms=7,
		label='active residual (+Ca, 200 s)')
	ax.axhspan(230, 270, color='#bdc3c7', alpha=0.35, label=r'target rest (~250 $\mu$M)')
	ax.axvline(0.135, color='#27ae60', ls='--', lw=1.3, label='v0.7 (0.135 pS)')
	ax.set_xlabel(r'$\gamma_{IP3R}$ (pS)')
	ax.set_ylabel(r'DTS store $\mathrm{Ca^{2+}}$ ($\mu$M)')
	ax.set_title('γ trades resting stability against active residual', fontsize=9.5)
	ax.legend(frameon=False, fontsize=8, loc='center right')
	ax.grid(alpha=0.3)

	fig.suptitle('Calibration: why V_IM = 0 mV and γ_IP3R = 0.135 pS', fontsize=12)
	fig.text(0.5, -0.05,
		'LEFT — sweeping V_IM (at the recalibrated γ) only moves the depletion '
		'floor: even at 0 mV the store still depletes ~99 % under sustained agonist, '
		'and the Dolan 5/5 peak/SOCE criteria pass at every V_IM (the cytosolic peak '
		'is buffer/PMCA/SOCE-clamped). V_IM therefore controls whether the store may '
		'fall below cytosol, not how deep it goes. RIGHT — sweeping γ (at V_IM = 0) '
		'shows the coupling: low γ under-leaks so the resting store over-fills '
		'(red rises far above the ~250 µM target); raising γ to ~0.135 pS lands the '
		'resting fixed point on target, but that same γ deepens active depletion, so '
		'the sustained-agonist residual (purple) collapses toward the cytosolic '
		'floor. A held residual is therefore unachievable by IP3R conductance alone '
		'— deep depletion under sustained saturating agonist is the honest outcome.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.05, 1, 0.95])
	path = os.path.join(out_dir, 'dts_calibration_sweeps.png')
	fig.savefig(path, dpi=150, bbox_inches='tight')
	plt.close(fig)
	print(f'  wrote {path}')


# ── Figure 4: dose-response depletion ─────────────────────────────────────

def fig_dose_response(out_dir: str) -> None:
	"""ADP dose-response (v0.7): the autocrine loop commits the store all-or-none,
	while the loops-off machinery shows a dose-graded residual."""
	doses = np.geomspace(0.05, 30.0, 9)
	full_resid, full_pk, lo_resid, lo_pk = [], [], [], []
	for dose in doses:
		f = run_with(*NEW, length=200, adp_peak_uM=float(dose), thrombin_off=True)
		l = run_with(*NEW, length=200, adp_peak_uM=float(dose), thrombin_off=True,
			loops_off=True)
		full_resid.append(f['ca_dts_uM'][-1]); full_pk.append(f['ca_cyt_nM'].max())
		lo_resid.append(l['ca_dts_uM'][-1]);   lo_pk.append(l['ca_cyt_nM'].max())

	fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
	ax = axes[0]
	ax.semilogx(doses, full_resid, 'o-', color='#c0392b', lw=2, ms=7,
		label='full model (autocrine loops on)')
	ax.semilogx(doses, lo_resid, 's--', color='#2980b9', lw=2, ms=6,
		label='loops off (Dolan-equivalent)')
	ax.set_xlabel(r'ADP dose ($\mu$M, log)')
	ax.set_ylabel(r'DTS residual at 200 s ($\mu$M)')
	ax.set_title('Sustained store residual vs dose', fontsize=10)
	ax.legend(frameon=False, fontsize=8, loc='upper right')
	ax.grid(alpha=0.3, which='both')

	ax = axes[1]
	ax.semilogx(doses, full_pk, 'o-', color='#c0392b', lw=2, ms=7,
		label='full model')
	ax.semilogx(doses, lo_pk, 's--', color='#2980b9', lw=2, ms=6, label='loops off')
	ax.axhline(200, color='#7f8c8d', ls=':', lw=1.2, label='Dolan "active" (200 nM)')
	ax.set_xlabel(r'ADP dose ($\mu$M, log)')
	ax.set_ylabel(r'peak cytosolic $\mathrm{Ca^{2+}}$ (nM)')
	ax.set_title(r'Peak cytosolic $\mathrm{Ca^{2+}}$ vs dose', fontsize=10)
	ax.legend(frameon=False, fontsize=8, loc='lower right')
	ax.grid(alpha=0.3, which='both')

	fig.suptitle('Dose-response: the autocrine loop is a commitment switch; the '
		'underlying store response is dose-graded', fontsize=12)
	fig.text(0.5, -0.05,
		'ADP-only dose-response of the v0.7 model (+Ca2+). In the full model (red) '
		'the autocrine ADP loop commits the store all-or-none: it empties and stays '
		'empty (residual ~0) at every supra-threshold dose, and the cytosolic peak '
		'is saturated across the range. Disable the autocrine loops (blue, '
		'Dolan-equivalent) and the underlying Ca-handling machinery is dose-graded — '
		'the store refills to a substantial residual at threshold (~44 µM at 0.05 µM '
		'ADP) and depletes toward empty as the dose saturates, and the peak grades '
		'down at the lowest doses. So deep depletion under a saturating agonist is '
		'the high-dose limit of a graded response; the autocrine loop just makes the '
		'platelet commit early. (The v0.7 γ recalibration deepened depletion overall, '
		'so the loops-off threshold residual is lower than in v0.63 — the store '
		'commitment is correspondingly sharper.)',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.05, 1, 0.95])
	path = os.path.join(out_dir, 'dose_response_depletion.png')
	fig.savefig(path, dpi=150, bbox_inches='tight')
	plt.close(fig)
	print(f'  wrote {path}')


def main(argv: list[str] | None = None) -> None:
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument('--out-dir', default='reports/figures/v0.7',
		help='Figure output dir. Default reports/figures/v0.7')
	args = p.parse_args(argv)
	fp.makedirs(args.out_dir)
	print(f'DTS-realism figures → {args.out_dir}/')
	fig_floor(args.out_dir)
	fig_resting(args.out_dir)
	fig_sweeps(args.out_dir)
	fig_dose_response(args.out_dir)
	print('done.')


if __name__ == '__main__':
	main()
