"""Unit tests for the CalciumTrace listener.

CalciumTrace reads from ``BulkMolecules.container.counts()`` and writes
columns of the Ca²⁺ signalling state for post-simulation analysis. These
tests bypass the simulation scaffold by setting the listener's indices
and bulk-molecules mock directly, so each test runs in <10 ms.

The SOCE flux calculation is independently re-implemented inside the
test to catch silent drift in the conversion factors (``_UM_PER_COUNT_*``)
or the Dolan Eq. 2–4 chain. The EDTA case (``CA_EX_UM = 0``) is exercised
to lock in the documented gating behaviour.
"""

import math
import unittest
from unittest.mock import MagicMock

import numpy as np

from models.platelet.listeners.calcium_trace import CalciumTrace
from reconstruction.platelet.dataclasses.process import calcium_signalling as cs_mod
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
	_UM_PER_COUNT_CYT,
	_UM_PER_COUNT_DTS,
	_mwc_open_fraction,
	GAMMA_SOC_S,
	NA_OVER_zF,
	ORAI_SUBUNITS_PER_CHANNEL,
	PUNCTA,
	RT_OVER_zF_V,
	V_PM_V,
)


# Order: ca_cyt, ca_dts, ip3, stim1, orai, cam_free, ca2_cam, ca4_cam,
#        ca4_cam_pmca, ca4_cam_pmca_ca, pmca_cam, pmca_free, pmca_ca
_NUM_SPECIES = 13


def _make_listener(counts):
	"""Build a CalciumTrace whose ``_bulk_molecules.counts()`` returns
	the supplied array, with the 13 species indices wired to slots 0–12.
	"""
	if len(counts) != _NUM_SPECIES:
		raise ValueError(
			f'counts must have length {_NUM_SPECIES}; got {len(counts)}')

	lst = CalciumTrace()
	bulk = MagicMock()
	bulk.counts.return_value = np.asarray(counts, dtype=np.int64)
	lst._bulk_molecules = bulk

	# CalciumDynamics reference, normally set in initialize(); the per-step
	# ATP cost is read from its `_atp_cost` attribute.
	cd = MagicMock()
	cd._atp_cost = 0
	lst._calcium_dynamics = cd

	lst._idx_ca_cyt          = 0
	lst._idx_ca_dts          = 1
	lst._idx_ip3             = 2
	lst._idx_stim1           = 3
	lst._idx_orai            = 4
	lst._idx_cam_free        = 5
	lst._idx_ca2_cam         = 6
	lst._idx_ca4_cam         = 7
	lst._idx_ca4_cam_pmca    = 8
	lst._idx_ca4_cam_pmca_ca = 9
	lst._idx_pmca_cam        = 10
	lst._idx_pmca_free       = 11
	lst._idx_pmca_ca         = 12
	return lst


class TestCalciumTraceUnitsConversion(unittest.TestCase):
	"""Direct pass-through fields apply the correct scaling factor."""

	def test_ca_cyt_nM_uses_cyt_factor_and_uM_to_nM(self):
		counts = [361] + [0] * (_NUM_SPECIES - 1)
		expected_nM = 361 * _UM_PER_COUNT_CYT * 1e3
		lst = _make_listener(counts)
		lst.update()
		self.assertAlmostEqual(lst.ca_cyt_nM, expected_nM, places=6)

	def test_ca_dts_uM_uses_dts_factor(self):
		counts = [0, 38_800] + [0] * (_NUM_SPECIES - 2)
		expected_uM = 38_800 * _UM_PER_COUNT_DTS
		lst = _make_listener(counts)
		lst.update()
		self.assertAlmostEqual(lst.ca_dts_uM, expected_uM, places=6)

	def test_ip3_nM_uses_cyt_factor(self):
		counts = [0, 0, 180] + [0] * (_NUM_SPECIES - 3)
		expected_nM = 180 * _UM_PER_COUNT_CYT * 1e3
		lst = _make_listener(counts)
		lst.update()
		self.assertAlmostEqual(lst.ip3_nM, expected_nM, places=6)


class TestCalciumTracePassThroughCounts(unittest.TestCase):
	"""CaM ladder + PMCA sub-state fields are integer copies of counts."""

	def test_all_cam_and_pmca_substate_fields(self):
		# ca_c, dts, ip3, stim, orai, cam_fr, ca2cam, ca4cam, cppm, cppca, pcam, pf, pca
		counts = [100, 38_800, 50, 200, 400,
			1_000, 500, 200,
			7, 8, 9, 30, 40]
		lst = _make_listener(counts)
		lst.update()

		self.assertEqual(lst.stim1_dim, 200)
		self.assertEqual(lst.cam_free, 1_000)
		self.assertEqual(lst.ca2_cam, 500)
		self.assertEqual(lst.ca4_cam, 200)
		self.assertEqual(lst.ca4_cam_pmca, 7)
		self.assertEqual(lst.ca4_cam_pmca_ca, 8)
		self.assertEqual(lst.pmca_cam, 9)
		self.assertEqual(lst.pmca_free, 30)
		self.assertEqual(lst.pmca_ca, 40)

	def test_atp_pump_per_s_reads_calcium_dynamics_cost(self):
		"""atp_pump_per_s mirrors CalciumDynamics._atp_cost as a float."""
		lst = _make_listener([0] * _NUM_SPECIES)
		lst._calcium_dynamics._atp_cost = 137
		lst.update()
		self.assertEqual(lst.atp_pump_per_s, 137.0)
		self.assertIsInstance(lst.atp_pump_per_s, float)


class TestCalciumTraceSoceFlux(unittest.TestCase):
	"""SOCE flux follows Dolan Eq. 2–4; zero under EDTA."""

	def setUp(self):
		self._saved_ca_ex_uM = cs_mod.CA_EX_UM

	def tearDown(self):
		cs_mod.CA_EX_UM = self._saved_ca_ex_uM

	def test_soce_flux_zero_when_ca_ex_zero(self):
		"""EDTA condition (``CA_EX_UM = 0``) ⇒ SOCE current = 0 nM/s."""
		cs_mod.CA_EX_UM = 0.0

		counts = [361, 38_800, 180, 200, 400] + [0] * (_NUM_SPECIES - 5)
		lst = _make_listener(counts)
		lst.update()

		self.assertEqual(lst.soce_flux_nMs, 0.0)

	def test_soce_flux_zero_when_ca_cyt_zero(self):
		"""ca_cyt = 0 short-circuits the Nernst E_Ca,PM (log(C_ex/0))."""
		cs_mod.CA_EX_UM = 1200.0
		counts = [0, 38_800, 180, 200, 400] + [0] * (_NUM_SPECIES - 5)
		lst = _make_listener(counts)
		lst.update()
		self.assertEqual(lst.soce_flux_nMs, 0.0)

	def test_soce_flux_matches_manual_hill_mwc_nernst_chain(self):
		"""Listener result agrees with a hand-replayed Dolan Eq. 2–4 chain."""
		cs_mod.CA_EX_UM = 1200.0
		ca_cyt_count = 361
		ca_dts_count = 38_800
		ip3_count = 180
		stim1_count = 200
		orai_count = 400
		counts = [ca_cyt_count, ca_dts_count, ip3_count,
			stim1_count, orai_count] + [0] * (_NUM_SPECIES - 5)

		lst = _make_listener(counts)
		lst.update()

		# Manual replay (kept independent of listener internals).
		ca_cyt_uM = ca_cyt_count * _UM_PER_COUNT_CYT
		hill = (ca_cyt_uM ** PUNCTA['n']
			/ (PUNCTA['KM_uM'] ** PUNCTA['n'] + ca_cyt_uM ** PUNCTA['n']))
		qp = PUNCTA['alpha'] * hill + PUNCTA['baseline']
		stim2_p = qp * stim1_count
		n_orai_channels = orai_count / ORAI_SUBUNITS_PER_CHANNEL
		po_orai, _ = _mwc_open_fraction(stim2_p, n_orai_channels)
		e_ca_pm_v = RT_OVER_zF_V * math.log(cs_mod.CA_EX_UM / ca_cyt_uM)
		driving_pm_v = V_PM_V - e_ca_pm_v
		soce_ions_s = (
			-GAMMA_SOC_S * n_orai_channels * po_orai
			* driving_pm_v * NA_OVER_zF)
		expected_nMs = soce_ions_s * _UM_PER_COUNT_CYT * 1e3

		self.assertAlmostEqual(lst.soce_flux_nMs, expected_nMs, places=8)

	def test_soce_flux_non_negative_at_resting_state(self):
		"""At resting cyt with extracellular > intracellular, SOCE ≥ 0."""
		cs_mod.CA_EX_UM = 1200.0
		# Resting-ish counts: cyt ~100 nM, basal STIM1, full Orai.
		counts = [361, 38_800, 180, 50, 400] + [0] * (_NUM_SPECIES - 5)
		lst = _make_listener(counts)
		lst.update()
		self.assertGreaterEqual(
			lst.soce_flux_nMs, 0.0,
			'SOCE should be a non-negative influx when extracellular '
			'Ca²⁺ exceeds cytosolic.')


if __name__ == '__main__':
	unittest.main()
