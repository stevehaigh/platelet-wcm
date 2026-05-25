"""Unit tests for the RestingDecay process.

RestingDecay implements stochastic exponential protein decay
(``lost ~ Binomial(count, p)``). These tests bypass the wholecell
simulation scaffold by instantiating the process directly and mocking
its bulk-molecules view, so each test runs in <10 ms.

The process uses ``np.random.binomial`` (NOT ``sim.randomState``), so
tests seed numpy's global RNG via ``np.random.seed``.
"""

import unittest
from unittest.mock import MagicMock

import numpy as np

from models.platelet.processes.resting_decay import RestingDecay


def _make_process(half_life_sec, dt_sec, counts):
	"""Build a RestingDecay with mocked timestep + bulk-molecules view.

	Returns (process, mock_view) so the caller can inspect calls to
	``countsDec``. The view's ``counts()`` returns the supplied array;
	mutate ``mock_view.counts.return_value`` between steps to model
	multi-step trajectories.
	"""
	proc = RestingDecay()
	proc._half_life = float(half_life_sec)

	view = MagicMock()
	view.counts.return_value = np.asarray(counts, dtype=np.int64)
	proc._molecules = view

	# Replace the bound timeStepSec() with a constant-returning callable.
	proc.timeStepSec = lambda: dt_sec  # type: ignore[assignment]
	return proc, view


class TestRestingDecayEdgeCases(unittest.TestCase):
	"""Trivial-answer cases — no RNG variance involved."""

	def test_zero_counts_yields_zero_loss(self):
		"""All-zero counts → countsDec called with all-zero loss vector."""
		proc, view = _make_process(
			half_life_sec=604800.0, dt_sec=1.0, counts=[0, 0, 0])

		np.random.seed(0)
		proc.evolveState()

		view.countsDec.assert_called_once()
		losses = view.countsDec.call_args[0][0]
		np.testing.assert_array_equal(losses, np.zeros(3, dtype=np.int64))

	def test_zero_dt_yields_zero_loss(self):
		"""dt = 0 → decay_prob = 0 → no loss regardless of count."""
		proc, view = _make_process(
			half_life_sec=604800.0, dt_sec=0.0,
			counts=[10_000, 50_000, 1])

		np.random.seed(123)
		proc.evolveState()

		losses = view.countsDec.call_args[0][0]
		np.testing.assert_array_equal(losses, np.zeros(3, dtype=np.int64))

	def test_infinite_half_life_yields_zero_loss(self):
		"""half_life → ∞ → decay_prob → 0 → no loss."""
		proc, view = _make_process(
			half_life_sec=np.inf, dt_sec=1.0, counts=[1_000_000])

		np.random.seed(42)
		proc.evolveState()

		losses = view.countsDec.call_args[0][0]
		np.testing.assert_array_equal(losses, np.zeros(1, dtype=np.int64))

	def test_countsDec_receives_int_array_of_same_shape(self):
		"""Loss vector must have the same shape as counts and integer dtype."""
		proc, view = _make_process(
			half_life_sec=604800.0, dt_sec=1.0,
			counts=[100, 200, 300, 400])

		np.random.seed(0)
		proc.evolveState()

		losses = view.countsDec.call_args[0][0]
		self.assertEqual(losses.shape, (4,))
		self.assertTrue(np.issubdtype(losses.dtype, np.integer))

	def test_calculate_request_calls_requestAll(self):
		"""calculateRequest must declare a request for all of its molecules."""
		proc, view = _make_process(
			half_life_sec=604800.0, dt_sec=1.0, counts=[1, 2, 3])

		proc.calculateRequest()

		view.requestAll.assert_called_once()


class TestRestingDecayBinomialStatistics(unittest.TestCase):
	"""The binomial draw must reproduce the analytical exponential on average."""

	def test_ensemble_mean_matches_expected_exponential(self):
		"""<N(t)> ≈ N₀ · exp(-ln2 · t / t_half) within 2 %.

		Ten 1-s timesteps with t_half = 50 s gives an expected survival
		fraction of ≈ 0.8706. With N₀ = 100 000 and 100 replicates the
		Monte-Carlo standard error is well under 0.5 %, so a ±2 % check
		is robust.
		"""
		half_life_sec = 50.0
		dt_sec = 1.0
		n_steps = 10
		n_replicates = 100
		n0 = 100_000

		expected_survival = np.exp(
			-np.log(2) * n_steps * dt_sec / half_life_sec)
		expected_count = n0 * expected_survival

		np.random.seed(7)
		final_counts = []
		for _ in range(n_replicates):
			proc, view = _make_process(half_life_sec, dt_sec, [n0])
			current = np.array([n0], dtype=np.int64)
			for _step in range(n_steps):
				view.counts.return_value = current
				proc.evolveState()
				losses = view.countsDec.call_args[0][0]
				current = current - losses
			final_counts.append(int(current[0]))

		mean_final = float(np.mean(final_counts))
		self.assertAlmostEqual(
			mean_final / expected_count, 1.0, delta=0.02,
			msg=(f'Ensemble mean {mean_final:.0f} deviates from expected '
				f'{expected_count:.0f} by >2 %.'))

	def test_half_life_equal_to_dt_gives_half_probability(self):
		"""dt = t_half → p = 1 - exp(-ln2) = 0.5; mean loss ≈ N₀/2."""
		half_life_sec = 1.0
		dt_sec = 1.0
		n0 = 1_000_000
		n_reps = 50

		np.random.seed(11)
		losses_each = []
		for _ in range(n_reps):
			proc, view = _make_process(half_life_sec, dt_sec, [n0])
			proc.evolveState()
			losses_each.append(int(view.countsDec.call_args[0][0][0]))

		mean_loss = float(np.mean(losses_each))
		self.assertAlmostEqual(
			mean_loss / (n0 * 0.5), 1.0, delta=0.01,
			msg=(f'Mean loss {mean_loss:.0f} ≠ N₀/2 ({n0 * 0.5:.0f}) '
				f'by >1 %.'))


if __name__ == '__main__':
	unittest.main()
