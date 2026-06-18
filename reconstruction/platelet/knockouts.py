"""
Logical-entity → species-id map for expression knockouts (Tier 2).

Zeroing a protein's resting copy number is an *expression knockout*. Several
proteins are tracked as multiple conformational sub-states (a receptor's
inactive / active / desensitised pool; SERCA's E1/E2 cycle), so a knockout must
zero the **whole** entity to preserve conservation — hence this map from one
logical name to all of its species ids. The ODE and processes read these
amounts from the state vector (e.g. ``v_par1_cleave = k_cleave · PAR1_inactive ·
[thrombin]``), so a count of 0 genuinely removes the protein's function.

Curated and extensible: add an entry only when a protein's dynamics are driven
by its state counts (mass-action), not by a TOML constant. The entries here are
all count-driven (verified against ``calcium_signalling._ode_rhs`` /
``IntegrinActivation``). ``RunConfig.count_overrides`` consumes the expanded
``{species_id: 0}`` mapping; the TUI offers these names as tick-box knockouts.
"""

from typing import Dict, Iterable, List

ENTITIES: Dict[str, List[str]] = {
	'P2Y1 (ADP receptor)': [
		'P2Y1_inactive[pl]', 'P2Y1_active[pl]', 'P2Y1_desensitised[pl]'],
	'PAR1 (thrombin)': [
		'PAR1_inactive[pl]', 'PAR1_active[pl]', 'PAR1_internalized[pl]'],
	'PAR4 (thrombin)': [
		'PAR4_inactive[pl]', 'PAR4_active[pl]', 'PAR4_internalized[pl]'],
	'P2X1 (ATP channel)': [
		'P2X1[pl]', 'P2X1_O[pl]', 'P2X1_D[pl]'],
	'TP (TXA2 receptor)': [
		'TP_inactive[pl]', 'TP_active[pl]'],
	'aIIbb3 (integrin)': [
		'aIIbb3_resting[pl]', 'aIIbb3_active[pl]'],
	'SERCA (DTS Ca2+ pump)': [
		'SERCA_E1[dts]', 'SERCA_E2[dts]', 'SERCA_E1Ca[dts]',
		'SERCA_E1PCa[dts]', 'SERCA_E2PCa[dts]', 'SERCA_E2P[dts]'],
}


def expand(entity_names: Iterable[str]) -> Dict[str, int]:
	"""Expand selected entity names into a ``{species_id: 0}`` knockout mapping."""
	overrides: Dict[str, int] = {}
	for name in entity_names:
		for species_id in ENTITIES[name]:
			overrides[species_id] = 0
	return overrides
