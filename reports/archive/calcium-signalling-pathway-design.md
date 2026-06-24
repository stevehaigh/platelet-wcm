# Modelling ADP Receptor Signalling to Ca2+ in the Platelet WCM

## Biology overview

**P2Y1 (Gq-coupled)** is the primary driver of Ca2+ release from intracellular stores. **P2Y12 (Gi-coupled)** amplifies and sustains the Ca2+ signal by inhibiting adenylyl cyclase, lowering cAMP, and reducing PKA-mediated inhibition of Ca2+ mobilisation. Both pathways are triggered by ADP and work in concert.

## Architecture pattern

The wcEcoli codebase already solves exactly this kind of problem. The **TwoComponentSystem** is the template -- it models a multi-step receptor -> kinase -> response regulator signalling cascade using:

1. **A stoichiometry matrix** defining all species and reactions
2. **Sympy-generated ODEs** with automatic Jacobian computation
3. **`scipy.integrate.solve_ivp`** (BDF/LSODA) to integrate within each timestep
4. **Integer molecule count <-> concentration conversion** at the boundary

The pattern splits across two files:

- **Dataclass** (`reconstruction/platelet/dataclasses/process/`) -- builds the stoich matrix, rate constants, ODE functions at parameter-fitting time
- **Process** (`models/platelet/processes/`) -- runs the ODE each timestep within the `calculateRequest()`/`evolveState()` lifecycle

## Species list for ADP -> Ca2+

Drawing from Purvis 2008 (upstream) and Dolan & Diamond 2014 (downstream):

### P2Y1 pathway (Gq -> Ca2+ release)

| Reaction | Description |
|----------|-------------|
| ADP[e] + P2Y1[m] <=> ADP:P2Y1[m] | Agonist binding |
| ADP:P2Y1[m] + Gaq-GDP[m] -> ADP:P2Y1:Gaq-GTP[m] | G-protein activation |
| Gaq-GTP[m] + PLCb[c] -> Gaq-GTP:PLCb*[c] | PLC activation |
| PIP2[m] --PLCb*--> IP3[c] + DAG[m] | Hydrolysis (catalytic) |
| IP3[c] + IP3R[dts] <=> IP3:IP3R*[dts] | IP3 receptor binding |
| Ca2+[dts] --IP3R*--> Ca2+[c] | Store release |
| Ca2+[c] --SERCA--> Ca2+[dts] | Reuptake |
| Ca2+[c] --PMCA--> Ca2+[e] | Extrusion |

### P2Y12 pathway (Gi -> amplification)

| Reaction | Description |
|----------|-------------|
| ADP[e] + P2Y12[m] <=> ADP:P2Y12[m] | Agonist binding |
| ADP:P2Y12[m] + Gai-GDP[m] -> Gai-GTP[m] + Gbg[m] | G-protein activation |
| Gai-GTP[m] + AC[m] -> Gai-GTP:AC[m] | AC inhibition |
| AC[m]: ATP[c] -> cAMP[c] | Reduced rate when inhibited |
| cAMP[c] + PKA[c] <=> cAMP:PKA*[c] | PKA activation |
| PKA*[c] --inhibits--> IP3R[dts] | PKA phosphorylates IP3R (reduces Ca2+ release) |

### SOCE (store-operated Ca2+ entry)

| Reaction | Description |
|----------|-------------|
| STIM1[dts] --senses low Ca2+[dts]--> STIM1*[m] | Oligomerisation |
| STIM1*[m] + Orai1[m] -> STIM1*:Orai1*[m] | CRAC channel opening |
| Ca2+[e] --Orai1*--> Ca2+[c] | Store-operated entry |

### Compartments

- `[c]` cytoplasm
- `[m]` membrane
- `[dts]` dense tubular system (ER equivalent)
- `[e]` extracellular

## Implementation in 3 layers

### Milestone 1 shortcut: hardcoded constants

For Milestone 1 we hardcode species and reaction data directly in the dataclass
module rather than reading TSV files.  This lets us validate the ODE biology
before investing in data pipeline boilerplate, and it **does not back us into a
corner** — the dataclass interface is identical either way.

The pattern is:

```python
# reconstruction/platelet/dataclasses/process/calcium_signalling.py

# ── Hardcoded data (replace with TSV parsing in Milestone 2) ──────────────
_SPECIES = [
    # (id,              compartment, initial_count, mass_fg)
    ('CA2[c]',          'c',         100,           0.0),
    ('CA2[dts]',        'dts',       11340,         0.0),
    ('IP3[c]',          'c',         0,             0.0),
    # ... full list from calcium-data-provenance.md
]

_REACTIONS = [
    # (reaction_id,         kf,      kr,    stoichiometry)
    ('IP3R_CA_RELEASE',  100.0, 0.0,
     {'IP3R_o[dts]': -1, 'CA2[dts]': -1, 'CA2[c]': +1}),
    ('SERCA_REUPTAKE',      4.0,     0.0,   {'CA2[c]': -1, 'CA2[dts]': +1}),
    # ... full list from calcium-data-provenance.md
]
# ──────────────────────────────────────────────────────────────────────────


class CalciumSignalling:
    def __init__(self, raw_data, sim_data):
        # Milestone 1: read hardcoded data above.
        # Milestone 2 migration: replace the two lines below with TSV parsing.
        #   species = _parse_species_tsv(raw_data.calcium_species)
        #   reactions = _parse_reactions_tsv(raw_data.calcium_reactions)
        species = _SPECIES
        reactions = _REACTIONS

        self._build_stoich_matrix(species, reactions)
        self._build_ode_system()
```

The process file (`models/platelet/processes/calcium_signalling.py`) and
`SimulationDataPlatelet` are **unchanged** by the migration.  The only edit
needed for Milestone 2 is swapping the two assignment lines above and adding the
TSV files.

### Milestone 2: raw data TSV files (deferred)

Create flat files analogous to `reconstruction/ecoli/flat/`:

```tsv
# reconstruction/platelet/flat/calcium_signalling_reactions.tsv
reaction_id       forward_rate  reverse_rate  stoichiometry
P2Y1_ADP_BINDING  1.0e7        10.0       ADP[e]:-1, P2Y1[m]:-1, ADP:P2Y1[m]:+1
GAQ_ACTIVATION    5.0e5         0.1           ADP:P2Y1[m]:-1, Gaq-GDP[m]:-1, ...
PLCB_ACTIVATION                1.0e6         1.0          ...
PIP2_HYDROLYSIS                50.0          0            ...
IP3R_BINDING                   8.0e6         2.0          ...
CA_RELEASE                     100.0         0            ...
SERCA_REUPTAKE                 4.0           0            ...
PMCA_EXTRUSION                 2.5           0            ...
```

Rate constants from Purvis 2008 (Table 1) and Dolan & Diamond 2014 (Table 2).

### Dataclass (parameter fitting) — Layer 2

```python
# reconstruction/platelet/dataclasses/process/calcium_signalling.py

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from wholecell.utils import build_ode


class CalciumSignalling:
    def __init__(self, raw_data, sim_data):
        # Parse stoichiometry from raw TSV
        # Build stoich matrix S (n_species x n_reactions)
        # Store rate constants kf, kr

        self._build_stoich_matrix(raw_data)
        self._build_ode_system()

    def _build_stoich_matrix(self, raw_data):
        """Parse reactions into sparse stoich matrix."""
        # Same pattern as TwoComponentSystem.__init__
        # molecules[], stoichMatrixI[], stoichMatrixJ[], stoichMatrixV[]
        ...

    def _build_ode_system(self):
        """Use sympy to build dy/dt = S . v(y) symbolically."""
        S = self.stoich_matrix()
        n_species, n_rxns = S.shape

        y = sp.symbols([f'y[{i}]' for i in range(n_species)])

        # Mass-action kinetics: v_j = kf * prod(reactants) - kr * prod(products)
        rates = []
        for j in range(n_rxns):
            fwd = self.rates_fwd[j]
            rev = self.rates_rev[j]
            for i in np.where(S[:, j] < 0)[0]:
                fwd *= y[i] ** (-S[i, j])
            for i in np.where(S[:, j] > 0)[0]:
                rev *= y[i] ** S[i, j]
            rates.append(fwd - rev)

        dy = sp·Matrix(S) * sp·Matrix(rates)
        J = dy.jacobian(y)

        # Compile to callable numpy functions
        self._derivatives, _ = build_ode.derivatives(dy)
        self._jacobian, _ = build_ode.derivatives_jacobian(J)

    def molecules_to_next_time_step(self, counts, volume,
                                     nAvogadro, dt, random_state):
        """Integrate ODEs for one timestep. Returns (needed, changes)."""
        y0 = counts / (volume * nAvogadro)  # counts -> concentration

        sol = solve_ivp(
            self._derivatives, [0, dt], y0,
            method='BDF', jac=self._jacobian,
            atol=1e-10, rtol=1e-8
        )

        y_final = np.maximum(sol.y[:, -1], 0)
        changes = np.round(
            (y_final - y0) * volume * nAvogadro
        ).astype(int)

        needed = (-changes).clip(min=0)
        return needed, changes
```

### Process (simulation runtime) — Layer 3

```python
# models/platelet/processes/calcium_signalling.py

import wholecell.processes.process as process
from wholecell.utils import units


class CalciumSignalling(process·Process):
    _name = "CalciumSignalling"

    def initialize(self, sim, sim_data):
        super().initialize(sim, sim_data)
        self.nAvogadro = sim_data.constants.n_avogadro.asNumber(1 / units.mol)
        self.cellDensity = sim_data.constants.cell_density.asNumber(
            units.g / units·L)

        # The ODE solver from the dataclass
        self.solver = sim_data.process.calcium_signalling

        # Build view over all species in the pathway
        self.molecules = self.bulkMoleculesView(
            sim_data.process.calcium_signalling.molecule_names)

    def calculateRequest(self):
        counts = self.molecules.total_counts()
        cellMass = (self.readFromListener("Mass", "cellMass")
                    * units.fg).asNumber(units.g)
        volume = cellMass / self.cellDensity

        self.needed, self.changes = (
            self.solver.molecules_to_next_time_step(
                counts, volume, self.nAvogadro,
                self.timeStepSec(), self.randomState))

        self.molecules.requestIs(self.needed)

    def evolveState(self):
        self.molecules.countsInc(self.changes)
```

The framework handles partitioning, mass tracking, and I/O around this.

## Key design decisions

### One process or many?

Model the entire Ca2+ signalling cascade as **one process**. The ODE system is tightly coupled (IP3 production drives Ca2+ release drives SERCA reuptake drives STIM1 sensing). Splitting it across processes would introduce artificial partitioning boundaries where the biology is continuous. The TwoComponentSystem does exactly this -- one process, one ODE system, many species.

### Where P2Y12 fits

P2Y12's Gi pathway modulates the P2Y1/Ca2+ pathway by reducing PKA inhibition. Include it as additional species and reactions in the same ODE system. The coupling is what makes the biology interesting -- P2Y12 knockout dampens the Ca2+ transient amplitude and duration, which the model should reproduce.

### Timestep considerations

Ca2+ transients in platelets peak in ~1-3 seconds and decay over ~30-60 seconds. With 1-second timesteps and BDF integration, this is well-resolved. If instability appears, the framework already supports `isTimeStepShortEnough()` to adaptively shorten.

### Validation target

The first milestone should be reproducing the Dolan & Diamond Ca2+ transient shape -- a sharp peak followed by a sustained plateau (SOCE-dependent). If the model matches their Figure 2 curves with the ODE parameters, that is a strong dissertation result on its own.

## Order of work

1. Start with the **Dolan & Diamond 34-species Ca2+ subsystem** alone (IP3R, SERCA, PMCA, SOCE) -- get the transient shape right
2. Add **Purvis upstream** (P2Y1 -> Gq -> PLCb -> IP3 generation) to make it agonist-driven
3. Add **P2Y12 -> Gi -> AC -> cAMP -> PKA** modulation to show the amplification effect
4. Validate: P2Y12 knockout (cangrelor) should reduce but not abolish the Ca2+ transient

Each step is independently testable and publishable.
