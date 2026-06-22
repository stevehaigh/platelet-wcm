"""
Shared helpers for the TUI per-theme demo figures (`demo_*.py`).

Each demo figure is focused on one theme (calcium / integrin / thromboxane /
secretion) so a given TUI demo shows only the panels relevant to it (see
`tui-demos.md`). They all support a grey **baseline overlay**: if the
environment variable ``PLATELET_BASELINE_SIMOUT`` points at another run's
``simOut`` directory, each trace is drawn over that run's trace in grey. The TUI
"Demo figure" button sets this to the run pinned via "Set baseline", so a
knockout / drug run reads directly against its intact control. Unset → no overlay.
"""

import os

from wholecell.io.tablereader import TableReader

BASELINE_ENV = 'PLATELET_BASELINE_SIMOUT'
BASELINE_COLOR = '0.7'  # grey


def read_col(sim_out_dir, listener, column):
	"""Read one listener column as a flat float array, or None if absent."""
	try:
		reader = TableReader(os.path.join(sim_out_dir, listener))
		return reader.readColumn(column).flatten().astype(float)
	except Exception:
		return None


def resolve_baseline(sim_out_dir):
	"""Baseline simOut dir from the env var, or None if unusable.

	None when unset, when it points at the current run (self-overlay is
	pointless), or when it holds no listener data.
	"""
	base = os.environ.get(BASELINE_ENV) or None
	if base is None:
		return None
	if os.path.abspath(base) == os.path.abspath(sim_out_dir):
		return None
	if not os.path.isdir(os.path.join(base, 'CalciumTrace')):
		return None
	return base


def make_draw(sim_out_dir, baseline_dir):
	"""Return a ``draw(ax, listener, column, color, label, ...)`` closure that
	plots the current run (coloured) over the baseline run (grey, behind)."""

	def draw(ax, listener, column, color, label, lw=2.0, ls='-', baseline=True):
		if baseline and baseline_dir is not None:
			bt = read_col(baseline_dir, listener, 'time')
			by = read_col(baseline_dir, listener, column)
			if bt is not None and by is not None and len(bt) == len(by) and len(bt):
				# Include the baseline run's duration so a length mismatch
				# (e.g. a 60 s baseline vs a longer current run) reads as
				# "the baseline run was shorter", not a broken/truncated line.
				blabel = ('%s (baseline · %.0f s)' % (label, bt[-1])
					if label else None)
				ax.plot(bt, by, color=BASELINE_COLOR, lw=1.3, ls=ls, zorder=1,
					label=blabel)
		t = read_col(sim_out_dir, listener, 'time')
		y = read_col(sim_out_dir, listener, column)
		if t is not None and y is not None and len(t) == len(y):
			ax.plot(t, y, color=color, lw=lw, ls=ls, zorder=2, label=label)

	return draw


def combined_legend(ax_left, ax_right, loc='best'):
	"""Merge the legends of a twin-axis pair onto the left axis."""
	h1, l1 = ax_left.get_legend_handles_labels()
	h2, l2 = ax_right.get_legend_handles_labels()
	ax_left.legend(h1 + h2, l1 + l2, loc=loc, fontsize=8)


def suptitle(fig, title, baseline_dir):
	"""Figure title, noting the grey baseline overlay when one is shown."""
	extra = '  ·  grey = baseline overlay' if baseline_dir is not None else ''
	fig.suptitle(title + extra, fontsize=12, fontweight='bold')
