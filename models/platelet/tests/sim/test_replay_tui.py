"""Smoke tests for the terminal replay TUI loader + rendering helpers.

The animation loop, keyboard listener, and rich Live display all need a
real TTY, so they're not exercised here. We do exercise:

- `_resolve_simout` (path-resolution edge cases)
- `load_snapshots` (full data load against a Phase-3-style sim)
- The Snapshot dataclass + Totals (every field populated, sane ranges)
- The pure-function rendering helpers (_bar, _sparkline, _ca_colour)

If `replayTui` is broken upstream of the rich rendering, these tests
will catch it.
"""

import os
import tempfile
import unittest

import numpy as np
import pytest

# The replayer is an optional visualisation extra (rich + textual,
# requirements-viz.txt). Skip the whole module cleanly if they're absent.
pytest.importorskip('rich', reason='TUI viz extra — see requirements-viz.txt')
pytest.importorskip('textual', reason='TUI viz extra — see requirements-viz.txt')

from runscripts.manual.replayTui import (
	Snapshot,
	Totals,
	_bar,
	_ca_colour,
	_resolve_simout,
	_sparkline,
	load_snapshots,
)
from runscripts.manual.runPlateletSim import run_platelet_sim


class TestReplayTuiHelpers(unittest.TestCase):
	"""Pure-function helpers — no IO."""

	def test_bar_full_and_empty(self):
		self.assertEqual(_bar(0.0, width=4), '▱▱▱▱')
		self.assertEqual(_bar(1.0, width=4), '▰▰▰▰')

	def test_bar_partial(self):
		# 0.5 of 4 → 2 full, 2 empty
		self.assertEqual(_bar(0.5, width=4), '▰▰▱▱')

	def test_bar_clamps(self):
		# Out-of-range inputs don't break.
		self.assertEqual(_bar(-1.0, width=4), '▱▱▱▱')
		self.assertEqual(_bar(2.0, width=4), '▰▰▰▰')

	def test_sparkline_constant_at_full_width(self):
		# Constant data renders as flat low blocks; must still be width-padded.
		spark = _sparkline(np.array([100.0] * 5), width=5)
		self.assertEqual(len(spark), 5)
		# Glyphs only, no padding (5 samples × width 5).
		self.assertNotIn(' ', spark)

	def test_sparkline_constant_short_data_pads(self):
		"""Constant series shorter than width must still be width-padded.

		Regression check for the panel-jitter bug Copilot caught on PR #48:
		early in a replay (when history < panel width) the constant-data
		branch returned `len(tail)` chars with no padding, causing the
		sparkline panel to visibly shrink and grow.
		"""
		spark = _sparkline(np.array([100.0, 100.0, 100.0]), width=10)
		self.assertEqual(len(spark), 10)
		self.assertTrue(spark.startswith(' '), 'expected left padding')
		self.assertNotEqual(spark[-1], ' ', 'expected glyph at right edge')

	def test_sparkline_varying_short_data_pads(self):
		"""Varying short data is right-aligned (left-padded with spaces) so
		the newest sample sits at the right edge."""
		spark = _sparkline(np.array([1.0, 2.0, 3.0]), width=10)
		self.assertEqual(len(spark), 10)
		self.assertTrue(spark.startswith(' '), 'expected left padding')
		self.assertNotEqual(spark[-1], ' ', 'expected glyph at right edge')

	def test_ca_colour_thresholds(self):
		self.assertEqual(_ca_colour(100), 'green')
		self.assertEqual(_ca_colour(200), 'yellow')
		self.assertEqual(_ca_colour(350), 'orange1')
		self.assertEqual(_ca_colour(500), 'red')


@pytest.mark.slow
class TestReplayTuiLoad(unittest.TestCase):
	"""End-to-end load against a fresh 5-second sim."""

	@classmethod
	def setUpClass(cls):
		cls._tmpdir = tempfile.TemporaryDirectory()
		cls.paths = run_platelet_sim(
			cls._tmpdir.name,
			length_sec=5, seed=0, log_to_shell=False)
		cls.simout = cls.paths['sim_out_dir']

	@classmethod
	def tearDownClass(cls):
		cls._tmpdir.cleanup()

	def test_resolve_simout_from_simout_path(self):
		# Direct path to simOut/. `Path.resolve()` follows symlinks
		# (e.g. macOS's /var → /private/var); compare canonical forms.
		resolved = _resolve_simout(self.simout)
		self.assertEqual(os.path.realpath(str(resolved)),
			os.path.realpath(self.simout))

	def test_resolve_simout_from_run_dir(self):
		# A higher-level run dir resolves into its inner simOut/.
		resolved = _resolve_simout(self._tmpdir.name)
		self.assertTrue(str(resolved).endswith('simOut'))
		self.assertTrue((resolved / 'CalciumTrace').is_dir())

	def test_resolve_simout_missing_raises(self):
		"""Pointing at a dir with no simOut/ should raise, not silently pick."""
		with tempfile.TemporaryDirectory() as empty:
			with self.assertRaises(FileNotFoundError):
				_resolve_simout(empty)

	def test_resolve_simout_ambiguous_raises(self):
		"""Multiple simOut/ candidates under a parent → error, not first-wins.

		Mirrors the Phase 3 driver layout where `out/phase3_xxx/` contains
		both `with_ca/` and `no_ca/` siblings, each with a simOut/.
		"""
		with tempfile.TemporaryDirectory() as parent:
			for name in ('with_ca', 'no_ca'):
				dest = os.path.join(parent, name, 'simOut')
				os.makedirs(os.path.join(dest, 'CalciumTrace'))
				os.makedirs(os.path.join(dest, 'BulkMolecules'))
			with self.assertRaises(FileNotFoundError) as ctx:
				_resolve_simout(parent)
			self.assertIn('Ambiguous', str(ctx.exception))

	def test_load_snapshots_length(self):
		# 5-s sim → 6 snapshots (t = 0, 1, 2, 3, 4, 5 inclusive).
		snaps, totals, meta = load_snapshots(_resolve_simout(self.simout))
		self.assertEqual(len(snaps), 6)
		self.assertIsInstance(snaps[0], Snapshot)
		self.assertIsInstance(totals, Totals)
		self.assertIsInstance(meta, dict)

	def test_snapshot_sane_at_t0(self):
		snaps, totals, _ = load_snapshots(_resolve_simout(self.simout))
		s0 = snaps[0]
		# Resting Ca²⁺ within physiological window.
		self.assertGreater(s0.ca_cyt_nM, 80.0)
		self.assertLess(s0.ca_cyt_nM, 120.0)
		# Resting IP3 around 50 nM.
		self.assertGreater(s0.ip3_nM, 30.0)
		self.assertLess(s0.ip3_nM, 80.0)
		# DTS in the 200-300 µM resting band.
		self.assertGreater(s0.ca_dts_uM, 200.0)
		self.assertLess(s0.ca_dts_uM, 300.0)

	def test_totals_positive(self):
		_, totals, _ = load_snapshots(_resolve_simout(self.simout))
		self.assertGreater(totals.par1, 0)
		self.assertGreater(totals.par4, 0)
		self.assertGreater(totals.p2y1, 0)
		self.assertGreater(totals.p2x1, 0)
		self.assertGreater(totals.gq, 0)
		self.assertGreater(totals.plcb, 0)
		self.assertGreater(totals.pmca, 0)
		self.assertGreater(totals.calr, 0)

	def test_mito_count_loaded(self):
		"""ca_mito_count populated from BulkMolecules; non-negative."""
		snaps, _, _ = load_snapshots(_resolve_simout(self.simout))
		self.assertGreaterEqual(snaps[0].ca_mito_count, 0)
		# Field is an int, not a float concentration — the engine
		# doesn't track mito volume separately.
		self.assertIsInstance(snaps[0].ca_mito_count, int)

	def test_cell_schematic_renders(self):
		"""_cell_schematic produces a Panel without raising; the membrane
		lines have the expected receptor tags + flux tags; DTS and Mito
		boxes both have matching top and bottom borders."""
		from runscripts.manual.replayTui import _cell_schematic
		from rich.console import Console
		from rich.panel import Panel
		snaps, totals, _ = load_snapshots(_resolve_simout(self.simout))
		# Mid-trace frame to exercise non-trivial receptor occupancy.
		mid = snaps[len(snaps) // 2]
		panel = _cell_schematic(mid, totals, ca_ex_uM=1200.0)
		# Returns a Panel — rich.Live drops bare Text from a Layout
		# region, so the schematic must be wrapped in a Panel even
		# though the cell already has its own ╔═╗ border.
		self.assertIsInstance(panel, Panel)
		# Render to plain text so we can grep for structural markers.
		console = Console(width=110, record=True, force_terminal=False,
			color_system=None)
		console.print(panel)
		body = console.export_text()
		# Membrane lines: receptor tags and flux tags present.
		self.assertIn('[ PAR1', body, 'PAR1 receptor missing from membrane')
		self.assertIn('[ PAR4', body, 'PAR4 receptor missing from membrane')
		self.assertIn('[ P2Y1', body, 'P2Y1 receptor missing from membrane')
		self.assertIn('[ P2X1', body, 'P2X1 receptor missing from membrane')
		self.assertIn('SOCE', body)
		self.assertIn('PMleak', body)
		# DTS and Mito boxes have both ┌ and └ borders (top + bottom).
		self.assertGreaterEqual(body.count('┌'), 2,
			'expected two ┌ corners (DTS top + Mito top)')
		self.assertGreaterEqual(body.count('└'), 2,
			'expected two └ corners (DTS bottom + Mito bottom)')

	def test_running_peak_and_auc_monotonic(self):
		"""peak_ca_so_far must be non-decreasing; AUC must be non-decreasing."""
		snaps, _, _ = load_snapshots(_resolve_simout(self.simout))
		peaks = [s.peak_ca_so_far for s in snaps]
		aucs = [s.auc_above_rest for s in snaps]
		for i in range(1, len(snaps)):
			self.assertGreaterEqual(peaks[i], peaks[i - 1],
				'running peak should never decrease')
			self.assertGreaterEqual(aucs[i], aucs[i - 1] - 1e-9,
				'AUC above rest is a cumulative integral; must be monotone')


@pytest.mark.slow
class TestTextualApp(unittest.TestCase):
	"""Construct-and-introspect smoke test for the Textual app.

	Textual apps can't run interactively in CI (no tty), but constructing
	one + reading its BINDINGS / CSS catches mistakes like misspelled key
	names, missing action handlers, or CSS syntax errors before they
	silently break the live experience.
	"""

	@classmethod
	def setUpClass(cls):
		cls._tmpdir = tempfile.TemporaryDirectory()
		from runscripts.manual.runPlateletSim import run_platelet_sim
		cls.paths = run_platelet_sim(
			cls._tmpdir.name, length_sec=5, seed=0, log_to_shell=False)

	@classmethod
	def tearDownClass(cls):
		cls._tmpdir.cleanup()

	def test_app_constructs_and_binds(self):
		from runscripts.manual.replayTui import (
			PlateletReplayApp, _resolve_simout, load_snapshots)
		snaps, totals, meta = load_snapshots(_resolve_simout(self.paths['sim_out_dir']))
		app = PlateletReplayApp(snaps, totals, meta,
			initial_speed=0.2, start_frame=0)

		# Every action referenced in BINDINGS must resolve to a real method.
		bound_actions = {b.action for b in app.BINDINGS}
		for action in bound_actions:
			# Textual conventions: action_X exists as a method, OR is one of
			# the built-in App actions (e.g. "quit"). We at least check that
			# either form is reachable.
			method = f'action_{action}'
			self.assertTrue(
				hasattr(app, method) or hasattr(app.__class__, f'action_{action}'),
				f'BINDING action {action!r} has no handler')

		# CSS must be a non-empty string (Textual will raise on parse if not).
		self.assertGreater(len(app.CSS), 0)
		# Snapshot count round-trips through the app.
		self.assertEqual(app._n, len(snaps))


if __name__ == '__main__':
	unittest.main()
