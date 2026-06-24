"""Why the PKC feedbacks barely move cytosolic Ca2+ under a saturating agonist.

The platelet Ca2+ response is store-limited: the DTS empties on essentially
every supra-threshold stimulus, so the cytosolic transient is pinned by the
SOCE/pump balance regardless of how the feedbacks tune the upstream Gq drive.
This script makes that point with three figures, then the companion
``runSecondWave.py`` shows the weak-transient-agonist regime where the feedbacks
*do* move cytosolic Ca2+ (the second wave).

  1. brake_effect_on_ca.png   v0.6 PKC brakes on vs off -> cytosolic Ca2+ ~unchanged
  2. why_brake_invisible.png  the store empties (~10 s) before the brake engages (~15-30 s)
  3. amplifiers_saturating.png v0.61 amplifiers on vs off (standard 300 s): Ca2+ flat,
                              IP3 small, P2Y1 desensitisation clear

Conditions use the supported override knobs (no monkeypatching), via
``RunConfig``:
  v0.5  = brakes off, amplifiers off   (k_des = k_plcb_phos = 0; gains 0)
  v0.6  = brakes on,  amplifiers off
  v0.61 = brakes on,  amplifiers on    (defaults)

Figure conventions: matplotlib **mathtext** for chemical formulae (raw unicode
super/subscripts drop glyphs), detailed legends, and a takeaway caption.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotStoreLimitedFeedbacks.py \\
        [--outdir out/figures] [--figure brake|why|amplifiers|all]
"""

from __future__ import annotations

import argparse
import os
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from runscripts.manual.runPlateletSim import run_platelet_sim
from reconstruction.platelet.run_config import RunConfig
from wholecell.io.tablereader import TableReader


def _run(brakes, amplifiers, length, **agonist):
	"""Run one sim with the v0.6 brakes and v0.61 amplifiers toggled via RunConfig.

	brakes:     P2Y1 desensitisation + PLCβ phosphorylation (k_des / k_plcb_phos
	            scales). amplifiers: autocrine ADP + thromboxane loops (autocrine
	            ADP gain / cox1_factor). Each is 1.0 (on) or 0.0 (off).
	"""
	run_config = RunConfig(
		ca_ex_mM=1.2,
		k_des_scale=1.0 if brakes else 0.0,
		k_plcb_phos_scale=1.0 if brakes else 0.0,
		autocrine_adp_gain=1.0 if amplifiers else 0.0,
		cox1_factor=1.0 if amplifiers else 0.0,
		**agonist)
	td = tempfile.mkdtemp()
	paths = run_platelet_sim(td, length_sec=length, seed=0,
		log_to_shell=False, run_config=run_config)
	return paths['sim_out_dir']


def _cal(sim_out, column):
	return TableReader(os.path.join(sim_out, 'CalciumTrace')
		).readColumn(column).flatten()


# ── 1. v0.6 brakes on vs off -> cytosolic Ca2+ ──────────────────────────────

def fig_brake_effect(outdir: str) -> str:
	"""v0.5 (brakes off) vs v0.6 (brakes on) cytosolic Ca2+, two stimuli.

	Standard (thrombin-dominated) and ADP-only (which unmasks the ADP-arm P2Y1
	brake). In both the two traces are superimposed: the brakes barely change the
	store-limited Ca2+ amplitude — their effect lives in IP3 / receptor state.
	"""
	conditions: list[tuple[str, dict[str, float]]] = [
		('Standard transient\n(thrombin + ADP + ATP)', {}),
		('ADP-only\n(thrombin off, ATP off)',
			dict(thrombin_peak_nM=0.0, atp_ex_peak_uM=0.0)),
	]
	fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True,
		gridspec_kw=dict(height_ratios=[3, 1]))
	for j, (title, kw) in enumerate(conditions):
		off = _cal(_run(False, False, 300, **kw), 'ca_cyt_nM')   # v0.5
		on = _cal(_run(True, False, 300, **kw), 'ca_cyt_nM')     # v0.6
		t = np.arange(len(off))
		ax, axd = axes[0, j], axes[1, j]
		ax.plot(t, off, lw=2.6, color='#888888',
			label='v0.5 — PKC brakes off')
		ax.plot(t, on, lw=1.5, ls='--', color='#c0392b',
			label='v0.6 — PKC brakes on (P2Y1 desens + PLCb phos)')
		ax.set_title(title, fontsize=10)
		ax.grid(alpha=0.3)
		if j == 0:
			ax.set_ylabel(r'cytosolic $\mathrm{Ca^{2+}}$ (nM)')
			ax.legend(frameon=False, fontsize=8, loc='lower right')
		axd.axhline(0, color='k', lw=0.6)
		axd.plot(t, on - off, lw=1.4, color='#2c3e50')
		axd.fill_between(t, on - off, 0, color='#2c3e50', alpha=0.15)
		axd.set_xlabel('time (s)')
		axd.grid(alpha=0.3)
		if j == 0:
			axd.set_ylabel('v0.6 − v0.5 (nM)')
	fig.suptitle(r'v0.6 PKC brakes barely move cytosolic $\mathrm{Ca^{2+}}$ '
		r'(store-limited; +$\mathrm{Ca^{2+}}$ 1.2 mM, 300 s)', fontsize=12)
	fig.text(0.5, -0.01,
		'The DTS store empties either way, so the Ca2+ amplitude is set by the '
		'SOCE/pump balance, not by the receptor drive. The brakes act hard '
		'upstream (IP3 down, P2Y1 ~60-75% desensitised) but almost none of that '
		'reaches free cytosolic Ca2+ — the difference panels are sub-nanomolar.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.01, 1, 0.96])
	out = os.path.join(outdir, 'brake_effect_on_ca.png')
	fig.savefig(out, dpi=150, bbox_inches='tight')
	plt.close(fig)
	return out


# ── 2. Why the brake is invisible: timing ───────────────────────────────────

def fig_why_invisible(outdir: str) -> str:
	"""Store depletion vs PKC-brake engagement timing (standard activation).

	The store is >98% empty by ~10 s, but the PKC brake (P2Y1 desensitisation)
	only ramps after ~15 s — by which time IP3R has already drained the store.
	The brake then cuts IP3 ~25%, but the plateau is SOCE-set (store-depletion
	gated), not IP3-set, so cytosolic Ca2+ doesn't follow.
	"""
	off = _run(False, False, 300)                      # v0.5 (brakes off)
	on = _run(True, False, 300)                        # v0.6 (brakes on)
	t = np.arange(len(_cal(on, 'ca_cyt_nM')))

	fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(8.5, 9), sharex=True)
	a1.plot(t, _cal(off, 'ca_cyt_nM'), lw=2.5, color='#888888',
		label='v0.5 — brakes off')
	a1.plot(t, _cal(on, 'ca_cyt_nM'), lw=1.5, ls='--', color='#c0392b',
		label='v0.6 — brakes on')
	a1.set_ylabel(r'cytosolic $\mathrm{Ca^{2+}}$ (nM)')
	a1.set_title('What we see: the two traces are essentially identical',
		fontsize=10)
	a1.axvspan(0, 12, color='#3498db', alpha=0.07)
	a1.legend(frameon=False, fontsize=9, loc='lower right')
	a1.grid(alpha=0.3)

	a2.plot(t, _cal(on, 'ca_dts_uM'), lw=2.2, color='#8e44ad',
		label=r'DTS store $\mathrm{Ca^{2+}}$')
	a2.axhline(0.02 * 250, color='k', lw=0.6, ls=':')
	a2.annotate('2% of resting store', (150, 0.02 * 250 + 8), fontsize=8)
	a2.set_ylabel(r'DTS store $\mathrm{Ca^{2+}}$ ($\mu$M)')
	a2.set_title(r'Why (1): the store is $>$98% empty by ~10 s — '
		'before the brake engages', fontsize=10)
	a2.axvspan(0, 12, color='#3498db', alpha=0.07)
	a2.grid(alpha=0.3)

	a3b = a3.twinx()
	a3.plot(t, _cal(off, 'ip3_nM'), lw=2.2, color='#888888',
		label=r'$\mathrm{IP_3}$ — brakes off')
	a3.plot(t, _cal(on, 'ip3_nM'), lw=2.2, color='#c0392b',
		label=r'$\mathrm{IP_3}$ — brakes on (braked down ~25%)')
	a3.set_ylabel(r'$\mathrm{IP_3}$ (nM)')
	a3b.plot(t, _cal(on, 'p2y1_desensitised_frac'), lw=1.8, ls='-.',
		color='#e67e22', label='P2Y1 desensitised fraction')
	a3b.set_ylabel('P2Y1 desensitised fraction', color='#e67e22')
	a3b.tick_params(axis='y', labelcolor='#e67e22')
	a3.set_xlabel('time (s)')
	a3.set_title(r'Why (2): the brake engages after ~15 s and cuts '
		r'$\mathrm{IP_3}$ — but only SOCE sets the plateau now', fontsize=10)
	a3.legend(frameon=False, fontsize=8, loc='upper left')
	a3b.legend(frameon=False, fontsize=8, loc='lower right')
	a3.axvspan(0, 12, color='#3498db', alpha=0.07)
	a3.grid(alpha=0.3)

	fig.suptitle('Why the brake is invisible: it arrives after the store is '
		'already drained, and aims at the wrong valve', fontsize=11)
	fig.tight_layout(rect=[0, 0, 1, 0.96])
	out = os.path.join(outdir, 'why_brake_invisible.png')
	fig.savefig(out, dpi=150, bbox_inches='tight')
	plt.close(fig)
	return out


# ── 3. v0.61 amplifiers under a saturating agonist ──────────────────────────

def fig_amplifiers_saturating(outdir: str) -> str:
	"""v0.61 amplifiers on vs off at the standard saturating agonist (300 s).

	Same store-limit story for the *amplifying* loops: cytosolic Ca2+ is
	unchanged, IP3 differs a little (~4%), and P2Y1 desensitisation clearly
	(~20%, autocrine-ADP driven). Compare runSecondWave.py for the weak-agonist
	regime where the amplifiers do move Ca2+.
	"""
	v06 = _run(True, False, 300)                        # brakes on, amplifiers off
	v061 = _run(True, True, 300)                        # full v0.61
	t = np.arange(len(_cal(v061, 'ca_cyt_nM')))
	panels = [
		('ca_cyt_nM', r'cytosolic $\mathrm{Ca^{2+}}$ (nM)',
			r'No real difference ($<$0.1%) — store-limited'),
		('ip3_nM', r'$\mathrm{IP_3}$ (nM)',
			r'Small difference (~4%): amplifiers add $\mathrm{G_q}$ drive'),
		('p2y1_desensitised_frac', 'P2Y1 desensitised fraction',
			'Clear difference (~20%): autocrine ADP drives P2Y1'),
	]
	fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
	for ax, (col, ylab, title) in zip(axes, panels):
		ax.plot(t, _cal(v06, col), lw=2.6, color='#888888',
			label='v0.6 — amplifiers off (brakes on)')
		ax.plot(t, _cal(v061, col), lw=1.6, ls='--', color='#c0392b',
			label='v0.61 — full (autocrine ADP + TXA$_2$ loops)')
		ax.set_xlabel('time (s)')
		ax.set_ylabel(ylab)
		ax.set_title(title, fontsize=10)
		ax.grid(alpha=0.3)
		ax.legend(frameon=False, fontsize=8, loc='best')
	fig.suptitle('Standard saturating agonist, 300 s — v0.61 autocrine '
		'amplifiers vs v0.6: the effect is in IP3 / receptor state, not Ca2+',
		fontsize=12)
	fig.text(0.5, -0.02,
		'Under a saturating agonist the store empties, so the amplifiers do not '
		'change the cytosolic Ca2+ trace; they show in the upstream readouts '
		'(IP3, P2Y1 desensitisation) and in the entirely new outputs '
		'(secretion, TXA2). runSecondWave.py shows the weak-agonist regime where '
		'they do move cytosolic Ca2+ (the second wave).',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.02, 1, 0.95])
	out = os.path.join(outdir, 'amplifiers_saturating.png')
	fig.savefig(out, dpi=150, bbox_inches='tight')
	plt.close(fig)
	return out


_FIGURES = {'brake': fig_brake_effect, 'why': fig_why_invisible,
	'amplifiers': fig_amplifiers_saturating}


def main(argv: list[str] | None = None) -> None:
	p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
	p.add_argument('--outdir', default='out/figures',
		help='Directory to write PNGs into (default out/figures).')
	p.add_argument('--figure', choices=list(_FIGURES) + ['all'], default='all',
		help='Which figure to generate (default all).')
	args = p.parse_args(argv)
	os.makedirs(args.outdir, exist_ok=True)
	keys = list(_FIGURES) if args.figure == 'all' else [args.figure]
	for key in keys:
		print(f'[{key}] wrote {_FIGURES[key](args.outdir)}')


if __name__ == '__main__':
	main()
