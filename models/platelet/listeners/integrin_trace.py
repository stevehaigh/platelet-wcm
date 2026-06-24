"""
IntegrinTrace listener for the platelet whole-cell model (v0.61 §3).

Records αIIbβ3 (GPIIb-IIIa) inside-out activation each timestep — the
high-affinity conformational fraction that flow cytometry reports via the
activation-specific antibody PAC-1, and the αIIbβ3-antagonist / Glanzmann
perturbation target.

Columns written:
  time             — simulation time (s)
  aIIbb3_active    — high-affinity (PAC-1⁺) αIIbβ3 count ([pl])
  aIIbb3_resting   — low-affinity αIIbβ3 count ([pl])
  active_frac      — active / (active + resting) — the per-cell PAC-1 readout
  integrin_gate    — PKC* × Ca²⁺ inside-out activation gate value (0–1)
  pka_brake        — PKA dis-inhibition factor on the inside-out rate (≥1;
                     1.0 at rest, >1 as P2Y12/Gi lowers cAMP — the clopidogrel
                     mechanism, issue #10)
  akt_active       — lumped PI3K/Akt activity (0–1); tracks P2Y12 occupancy,
                     dis-inhibits the Rap1b-GAP (#73)
  rap1b_gtp        — Rap1b-GTP activity (0–1) — the integrin's proximal
                     forward driver; sustained by Akt while P2Y12 is driven (#73)
"""

import wholecell.listeners.listener


class IntegrinTrace(wholecell.listeners.listener.Listener):
	"""Record αIIbβ3 inside-out activation state each timestep."""

	_name = 'IntegrinTrace'

	def __init__(self, *args, **kwargs):
		self._bulk_molecules = None
		super().__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super().initialize(sim, sim_data)

		self._bulk_molecules = sim.internal_states['BulkMolecules'].container
		self._activation = sim.processes['IntegrinActivation']

		all_ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
		self._idx_active = all_ids.index('aIIbb3_active[pl]')
		self._idx_resting = all_ids.index('aIIbb3_resting[pl]')

		self.aIIbb3_active = 0
		self.aIIbb3_resting = 0
		self.active_frac = 0.0
		self.integrin_gate = 0.0
		self.pka_brake = 1.0
		self.akt_active = 0.0
		self.rap1b_gtp = 0.0

		self.registerLoggedQuantity('αIIbβ3 act\n(frac)', 'active_frac', '.3f')
		self.registerLoggedQuantity('PKA brake\n(≥1)', 'pka_brake', '.3f')

	def update(self):
		counts = self._bulk_molecules.counts()
		active = int(counts[self._idx_active])
		resting = int(counts[self._idx_resting])
		self.aIIbb3_active = active
		self.aIIbb3_resting = resting
		total = active + resting
		self.active_frac = float(active / total) if total > 0 else 0.0
		self.integrin_gate = float(getattr(self._activation, '_gate', 0.0))
		self.pka_brake = float(getattr(self._activation, '_pka_brake', 1.0))
		self.akt_active = float(getattr(self._activation, '_akt', 0.0))
		self.rap1b_gtp = float(getattr(self._activation, '_rap', 0.0))

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			units=('count for active/resting; active_frac, gate, pka_brake, '
				'akt_active, rap1b_gtp dimensionless'),
		)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time=self.time(),
			simulationStep=self.simulationStep(),
			aIIbb3_active=self.aIIbb3_active,
			aIIbb3_resting=self.aIIbb3_resting,
			active_frac=self.active_frac,
			integrin_gate=self.integrin_gate,
			pka_brake=self.pka_brake,
			akt_active=self.akt_active,
			rap1b_gtp=self.rap1b_gtp,
		)
