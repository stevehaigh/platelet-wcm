"""Figures for the v0.61 downstream PKC effects — what each new module does.

Three self-contained figures, one per v0.61 module. Each runs its own
simulation(s) (using the supported override knobs — no monkeypatching) and
writes a PNG to the output directory (default ``out/figures``).

  1. granule secretion       -> secretion_release.png
  2. autocrine ADP loop      -> autocrine_adp_loop.png
  3. thromboxane TXA2->TP->Gq -> thromboxane_loop.png

These show the new mechanisms *functioning*; the companion script
``plotStoreLimitedFeedbacks.py`` shows why they barely move cytosolic Ca2+
under a saturating agonist, and ``runSecondWave.py`` shows the weak-agonist
regime where they do.

Figure conventions (see the project figure-style notes): matplotlib **mathtext**
for chemical formulae (raw unicode super/subscripts drop glyphs in the default
font), detailed per-series legends, and a caption stating the takeaway.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotDownstreamModules.py \\
        [--outdir out/figures] [--figure secretion|autocrine|thromboxane|all]
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
import reconstruction.platelet.dataclasses.process.thromboxane_synthesis as tx_mod
from wholecell.io.tablereader import TableReader


def _run(length, **agonist):
	"""Run one platelet sim into a temp dir; return its simOut path."""
	td = tempfile.mkdtemp()
	paths = run_platelet_sim(td, length_sec=length, seed=0, log_to_shell=False,
		ca_ex_mM=1.2, **agonist)
	return paths['sim_out_dir']


def _col(sim_out, listener, column):
	return TableReader(os.path.join(sim_out, listener)).readColumn(column).flatten()


# ── 1. Granule secretion ────────────────────────────────────────────────────

def fig_secretion(outdir: str) -> str:
	"""Dense- vs alpha-granule cargo release under a standard activation.

	Shows: release is zero until PKC activates (~10 s), graded thereafter, and
	dense granules (ADP / serotonin) release faster than alpha granules
	(fibrinogen / P-selectin). The secretion gate (PKC x Ca2+) is overlaid.
	"""
	sim = _run(200)                                    # default agonist transient
	t = np.arange(len(_col(sim, 'SecretionTrace', 'secretion_gate')))
	adp = 100 * _col(sim, 'SecretionTrace', 'adp_released_frac')
	sht = 100 * _col(sim, 'SecretionTrace', 'serotonin_released_frac')
	fga = 100 * _col(sim, 'SecretionTrace', 'fibrinogen_released_frac')
	psel = 100 * _col(sim, 'SecretionTrace', 'pselectin_surface_frac')
	gate = _col(sim, 'SecretionTrace', 'secretion_gate')

	fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.7))
	a1.plot(t, adp, lw=2.4, color='#c0392b', label='ADP (dense granule)')
	a1.plot(t, sht, lw=1.4, ls='--', color='#e67e22',
		label='serotonin / 5-HT (dense granule)')
	a1.plot(t, fga, lw=2.4, color='#2980b9', label='fibrinogen (alpha granule)')
	a1.plot(t, psel, lw=1.4, ls='--', color='#27ae60',
		label='P-selectin surface (alpha granule)')
	a1.axhline(0, color='#999', lw=0.8, ls=':')
	a1.set_xlabel('time (s)')
	a1.set_ylabel('cargo released / surface-exposed (%)')
	a1.set_title('Granule cargo release (standard activation)', fontsize=10)
	a1.legend(frameon=False, fontsize=8, loc='center right')
	a1.grid(alpha=0.3)

	a2.plot(t, gate, lw=2.4, color='#8e44ad',
		label=r'secretion gate = PKC$^*$ $\times$ Ca$^{2+}$')
	a2.fill_between(t, gate, 0, color='#8e44ad', alpha=0.12)
	a2.set_xlabel('time (s)')
	a2.set_ylabel('secretion gate (0–1)')
	a2.set_title(r'Gate: 0 at rest (PKC below floor) $\rightarrow$ ramps on activation',
		fontsize=10)
	a2.legend(frameon=False, fontsize=9, loc='lower right')
	a2.grid(alpha=0.3)

	fig.suptitle(r'v0.61 granule secretion (+$\mathrm{Ca^{2+}}$ 1.2 mM, '
		'standard thrombin + ADP + ATP transient)', fontsize=12)
	fig.text(0.5, -0.02,
		'PKC + Ca2+-gated SNARE secretion relocates pre-existing cargo to the '
		'extracellular space (P-selectin to the surface). Release is exactly '
		'zero at rest (the gate keys off PKC activation above a resting-tone '
		'floor); dense-granule release leads, alpha-granule release follows.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.02, 1, 0.95])
	out = os.path.join(outdir, 'secretion_release.png')
	fig.savefig(out, dpi=150, bbox_inches='tight')
	plt.close(fig)
	return out


# ── 2. Autocrine ADP loop ───────────────────────────────────────────────────

def fig_autocrine(outdir: str) -> str:
	"""Thrombin-only run: secreted ADP closes the loop onto P2Y1.

	Thrombin (no exogenous ADP) drives PKC -> secretion. The secreted ADP is the
	*only* P2Y1 ligand, so any P2Y1 engagement here is purely autocrine. Shows
	the pericellular ADP rising then being cleared by ecto-NTPDase (AMP rises),
	and the P2Y1 it drives — a self-limiting loop.
	"""
	sim = _run(200, thrombin_peak_nM=1.0, adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
	t = np.arange(len(_col(sim, 'SecretionTrace', 'adp_e_uM')))
	adp_e = _col(sim, 'SecretionTrace', 'adp_e_uM')
	amp = _col(sim, 'SecretionTrace', 'amp_e')
	des = _col(sim, 'CalciumTrace', 'p2y1_desensitised_frac')

	fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.7))
	a1.plot(t, adp_e, lw=2.4, color='#c0392b',
		label=r'secreted ADP$_{[e]}$ (autocrine ligand)')
	a1.axhline(1.0, ls=':', color='#999', lw=1.1)
	a1.text(120, 1.05, r'$\approx$ P2Y1 threshold', color='#777', fontsize=8)
	a1.set_xlabel('time (s)')
	a1.set_ylabel(r'pericellular [ADP] ($\mu$M)', color='#c0392b')
	a1.tick_params(axis='y', labelcolor='#c0392b')
	a1b = a1.twinx()
	a1b.plot(t, amp / 1e3, lw=1.9, ls='--', color='#16a085',
		label=r'AMP$_{[e]}$ (ecto-NTPDase product)')
	a1b.set_ylabel(r'AMP$_{[e]}$ ($\times 10^3$ molecules)', color='#16a085')
	a1b.tick_params(axis='y', labelcolor='#16a085')
	a1.set_title('Secreted ADP rises, then is cleared (self-limiting)', fontsize=10)
	# Combined legend.
	h1, l1 = a1.get_legend_handles_labels()
	h2, l2 = a1b.get_legend_handles_labels()
	a1.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc='upper right')
	a1.grid(alpha=0.3)

	a2.plot(t, des, lw=2.4, color='#8e44ad',
		label='P2Y1 desensitised fraction (autocrine-driven)')
	a2.fill_between(t, des, 0, color='#8e44ad', alpha=0.12)
	a2.set_xlabel('time (s)')
	a2.set_ylabel('P2Y1 desensitised fraction')
	a2.set_title('P2Y1 driven entirely by autocrine ADP\n(thrombin-only: zero exogenous ADP)',
		fontsize=10)
	a2.legend(frameon=False, fontsize=8, loc='lower right')
	a2.grid(alpha=0.3)

	fig.suptitle(r'v0.61 autocrine ADP loop: thrombin $\rightarrow$ PKC '
		r'$\rightarrow$ secretion $\rightarrow$ ADP$_{[e]}$ $\rightarrow$ P2Y1',
		fontsize=12)
	fig.text(0.5, -0.02,
		'With no exogenous ADP, PKC-triggered dense-granule secretion is the only '
		'P2Y1 ligand. Secreted ADP reaches the P2Y1 range (~2 uM) and drives the '
		'receptor (engagement -> PKC-mediated desensitisation), then ecto-NTPDase '
		'hydrolyses it to AMP — so the loop self-limits even without exhausting '
		'the granule pool.',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.02, 1, 0.95])
	out = os.path.join(outdir, 'autocrine_adp_loop.png')
	fig.savefig(out, dpi=150, bbox_inches='tight')
	plt.close(fig)
	return out


# ── 3. Thromboxane TXA2 -> TP -> Gq ─────────────────────────────────────────

def fig_thromboxane(outdir: str) -> str:
	"""TXA2 synthesis + the autocrine TP->Gq loop, with the aspirin knockout.

	Standard activation, COX-1 intact vs aspirin (COX1_FACTOR = 0). Shows TXA2
	activating TP (left), and the loop's modest amplification of IP3 vs aspirin
	(right) — modest because the response is store-limited and PAR-dominated.
	"""
	cox0 = tx_mod.COX1_FACTOR
	try:
		tx_mod.COX1_FACTOR = 1.0
		on = _run(200)
		txa2 = _col(on, 'ThromboxaneTrace', 'txa2_uM')
		tp = _col(on, 'ThromboxaneTrace', 'tp_active_frac')
		ip3_on = _col(on, 'CalciumTrace', 'ip3_nM')
		tx_mod.COX1_FACTOR = 0.0
		asp = _run(200)
		ip3_asp = _col(asp, 'CalciumTrace', 'ip3_nM')
	finally:
		tx_mod.COX1_FACTOR = cox0
	t = np.arange(len(txa2))

	fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.7))
	a1.plot(t, txa2, lw=2.4, color='#c0392b', label=r'TXA$_2$ ($\mu$M)')
	a1.plot(t, tp, lw=2.4, color='#8e44ad', label='TP receptor active fraction')
	a1.set_xlabel('time (s)')
	a1.set_ylabel(r'TXA$_2$ ($\mu$M)  /  TP active fraction')
	a1.set_title(r'Synthesised TXA$_2$ activates the TP receptor', fontsize=10)
	a1.legend(frameon=False, fontsize=9, loc='lower right')
	a1.grid(alpha=0.3)

	a2.plot(t, ip3_on, lw=2.4, color='#16a085',
		label='COX-1 intact (TXA$_2$ loop ON)')
	a2.plot(t, ip3_asp, lw=1.9, ls='--', color='#7f8c8d',
		label='aspirin: COX-1 knockout (loop OFF)')
	a2.set_xlabel('time (s)')
	a2.set_ylabel(r'$\mathrm{IP_3}$ (nM)')
	a2.set_title(r'Autocrine amplification of the $\mathrm{G_q}$ cascade',
		fontsize=10)
	a2.legend(frameon=False, fontsize=9, loc='lower right')
	a2.grid(alpha=0.3)

	a3.plot(t, ip3_on - ip3_asp, lw=2.4, color='#16a085')
	a3.fill_between(t, ip3_on - ip3_asp, 0, color='#16a085', alpha=0.12)
	a3.axhline(0, color='#999', lw=0.8, ls=':')
	a3.set_xlabel('time (s)')
	a3.set_ylabel(r'$\Delta\mathrm{IP_3}$ (nM)')
	a3.set_title(r'$\Delta\mathrm{IP_3}$ from the loop (intact $-$ aspirin)',
		fontsize=10)
	a3.grid(alpha=0.3)

	fig.suptitle(r'v0.61 thromboxane loop: TXA$_2$ $\rightarrow$ TP '
		r'$\rightarrow$ G$_q$ (+ aspirin knockout)', fontsize=12)
	fig.text(0.5, -0.02,
		'TXA2 reaches ~1 uM and activates TP, which joins the shared Gq drive — '
		'closing the second autocrine amplifier. The effect on IP3 is modest '
		'(store-limited response, PAR-dominated under strong thrombin) and on '
		'cytosolic Ca2+ negligible; the third panel isolates the loop\'s IP3 '
		'contribution as the difference vs aspirin. Aspirin (COX-1 = 0) removes '
		'the loop entirely. The amplifier bites harder at a weak agonist (see '
		'runSecondWave.py).',
		ha='center', va='top', fontsize=8, wrap=True)
	fig.tight_layout(rect=[0, 0.02, 1, 0.95])
	out = os.path.join(outdir, 'thromboxane_loop.png')
	fig.savefig(out, dpi=150, bbox_inches='tight')
	plt.close(fig)
	return out


_FIGURES = {'secretion': fig_secretion, 'autocrine': fig_autocrine,
	'thromboxane': fig_thromboxane}


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
