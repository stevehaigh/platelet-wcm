# wcEcoli Reading Notes: How the Whole-Cell Model Works

A deep dive into the wcEcoli codebase—what it does, how it works, and what you can use it for.

**Key papers:**
- Macklin et al. (2020) "Simultaneous cross-evaluation of heterogeneous E. coli datasets via mechanistic simulation" *Science* 369:eaav3751
- Karr et al. (2012) "A Whole-Cell Computational Model Predicts Phenotype from Genotype" *Cell* 150(2):389-401

---

## The Shortest Mental Model

`wcEcoli` has three layers:

1. **Framework** in `wholecell/` — generic, organism-agnostic simulation engine
2. **E. coli model** in `models/ecoli/` — biological processes, listeners, analysis
3. **Parameter reconstruction** in `reconstruction/ecoli/` — raw data → fitted `sim_data`

The overall workflow:

1. Read experimental inputs and fit them into a simulation-data object
2. Run a simulation timestep by timestep using that fitted data
3. Record outputs to disk through listeners
4. Read those outputs back for analysis plots or web inspection

If you keep those four steps in mind, most of the repository layout makes sense.

---

## Start Here: First Files to Open

For a newcomer orienting themselves in the repository:

- `wholecell/sim/simulation.py` — the generic simulation engine and core loop
- `models/ecoli/sim/simulation.py` — E. coli-specific assembly: which states,
  processes, and listeners are active, and process dependency ordering
- `wholecell/processes/process.py` — base class for a process; the interface
  all biological submodels implement
- `reconstruction/ecoli/fit_sim_data_1.py` — the parameter calculator (ParCa);
  where raw knowledge-base data becomes the `sim_data` object
- `runscripts/manual/runParca.py` and `runscripts/manual/runSim.py` — main
  command-line entry points
- `wholecell/io/tablewriter.py` and `wholecell/io/tablereader.py` — the custom
  columnar output format; bridges simulation to analysis

The most useful habit in this repository is asking which layer a change
belongs to: framework, E. coli model, reconstruction, or analysis/UI. Many
confusing bugs become simpler once placed in the right layer.

---

## What Does wcEcoli Actually Do?

wcEcoli simulates a single *E. coli* cell from birth to division (approximately 60-120 minutes of biological time). It tracks:

- **Every gene** (~4,600 genes) and whether it's being transcribed
- **Every mRNA molecule** (~2,000-8,000 copies) and its degradation
- **Every protein** (~2.4 million total protein molecules) being synthesized and degraded
- **Every metabolite** (~1,500 species) concentrations and fluxes
- **DNA replication** fork positions and chromosome copy number
- **Cell growth** from ~650 fg dry mass to ~1,300 fg at division

The simulation predicts emergent phenotypes from genotype: growth rate, metabolic fluxes, gene expression levels, and cell cycle timing all emerge from the mechanistic models rather than being imposed.

---

## Inputs: What Data Goes In

### Raw Data Files (72+ TSV files in `reconstruction/ecoli/flat/`)

| Category | Example Files | What They Contain |
|----------|--------------|-------------------|
| **Genome** | `genes.tsv`, `sequence.fasta` | 4,583 genes, coordinates, strand, essential annotations |
| **Transcription** | `rnas.tsv`, `transcription_units.tsv`, `transcription_factors.tsv` | RNA properties, operons, TF binding sites |
| **Translation** | `proteins.tsv`, `translation_efficiency.tsv` | Protein sequences, molecular weights, translation rates |
| **Metabolism** | `metabolic_reactions.tsv`, `metabolites.tsv`, `metabolism_kinetics.tsv` | ~2,500 reactions, stoichiometry, kinetic parameters |
| **Regulation** | `fold_changes.tsv`, `ppgpp_regulation.tsv` | TF-gene fold changes, ppGpp effects |
| **Biomass** | `biomass.tsv`, `dry_mass_composition.tsv` | Target cell composition at different growth rates |

**Data format example** (`genes.tsv`):
```tsv
# Source: EcoCyc database, PMID: 12345678
id	name	symbol	left_end_pos	right_end_pos	strand	is_essential
EG10001	thrA	thrA	337	2799	+	True
```

### Parameter Calculator (ParCa) Processing

The raw data is transformed by `reconstruction/ecoli/fit_sim_data_1.py` through these stages:

```
Raw TSV files
  → initialize()           Load and validate raw data
  → input_adjustments()    Apply manual corrections
  → basal_specs()          Calculate basal expression levels
  → tf_condition_specs()   Fit TF regulation parameters  
  → fit_condition()        CVXPY optimization of expression
  → promoter_binding()     Calculate TF-promoter binding probabilities
  → final_adjustments()    Calibrate to target growth rate
  → sim_data pickle        ~180 MB serialized object
```

**Key fitted parameters:**
- RNA synthesis probabilities (which genes get transcribed and how often)
- Ribosome and RNAP expression levels (to achieve target growth rate)
- Metabolic enzyme kcat values (from literature + fitting)
- Translation efficiencies (from ribosome profiling data)

**ParCa runtime:** ~18 minutes, produces `out/{sim_dir}/kb/simData.cPickle`

---

## The Simulation Loop: What Happens Each Timestep

### Overview

Each simulation timestep (default 1 second, dynamically adjusted 0.1-1.0s):

```
┌─────────────────────────────────────────────────────────────┐
│                    TIMESTEP (1 second)                      │
├─────────────────────────────────────────────────────────────┤
│  1. CALCULATE REQUESTS                                      │
│     Each process declares what molecules it needs           │
│                                                             │
│  2. PARTITION STATE                                         │
│     Molecules allocated to processes by priority            │
│                                                             │
│  3. EVOLVE STATE                                            │
│     Processes execute their biological function             │
│                                                             │
│  4. MERGE STATE                                             │
│     Updates combined back into global state                 │
│                                                             │
│  5. UPDATE LISTENERS                                        │
│     Record data to disk for analysis                        │
└─────────────────────────────────────────────────────────────┘
```

### Process Execution Order

Processes run in **dependency-ordered groups** (from `models/ecoli/sim/simulation.py`):

```python
_processClasses = (
    # Group 1: TF unbinding (must run first)
    (TfUnbinding,),
    
    # Group 2: Chemical equilibria
    (Equilibrium, TwoComponentSystem, RnaMaturation,),
    
    # Group 3: TF binding (before transcription)
    (TfBinding,),
    
    # Group 4: Initiation + degradation (parallel)
    (TranscriptInitiation, PolypeptideInitiation, 
     ChromosomeReplication, ProteinDegradation, 
     RnaDegradation, Complexation,),
    
    # Group 5: Elongation (parallel)
    (TranscriptElongation, PolypeptideElongation,),
    
    # Group 6: Chromosome structure
    (ChromosomeStructure,),
    
    # Group 7: Metabolism (FBA, runs last to balance)
    (Metabolism,),
    
    # Group 8: Cell division check
    (CellDivision,),
)
```

**Why this order matters:**
- TF binding state must be resolved before transcription initiation can determine which genes are expressed
- Elongation must follow initiation
- Metabolism runs last to balance all the molecular changes from other processes
- Division checks if chromosome replication is complete

### State Partitioning: The Resource Competition Model

Processes compete for molecules through a priority-based allocation system:

```python
# Priority levels (lower = higher priority)
REQUEST_PRIORITY_METABOLISM = -10      # Gets first access to metabolites
REQUEST_PRIORITY_DEFAULT = 0           # Most processes
REQUEST_PRIORITY_DEGRADATION = 10      # Gets leftovers
```

**Partitioning algorithm:**
1. Each process requests N molecules of species X
2. Available molecules = total count - reserved for higher priority
3. If requests > available: proportional allocation based on request size
4. Processes receive their allocated share for use in `evolveState()`

---

## How Each Major Process Works

### Transcription (TranscriptInitiation + TranscriptElongation)

**TranscriptInitiation** (`models/ecoli/processes/transcript_initiation.py`):
```
Inputs:
  - Free RNA polymerase count
  - Promoter binding probabilities (from TF state)
  - Gene copy number (from replication state)
  - ppGpp concentration (stress signal)

Algorithm:
  1. Calculate synthesis probability for each gene:
     P(synthesis) = basal_prob × TF_modulation × ppGpp_effect
  2. Multinomial draw: which genes get initiated this timestep
  3. Create new "active_RNAP" unique molecules at gene start positions

Outputs:
  - Active RNAP molecules on DNA
  - Reduced free RNAP count
```

**TranscriptElongation**:
```
Inputs:
  - Active RNAP positions on DNA
  - NTP (nucleotide) availability
  - Elongation rate (~50 nt/s, varies with ppGpp)

Algorithm:
  1. For each RNAP: advance position by elongation_rate × timestep
  2. Consume NTPs proportional to sequence
  3. If reached gene end: create new RNA molecule, free RNAP

Outputs:
  - Updated RNAP positions
  - New RNA molecules (mRNA, rRNA, tRNA)
  - NTP consumption
```

### Translation (PolypeptideInitiation + PolypeptideElongation)

**PolypeptideElongation** (`models/ecoli/processes/polypeptide_elongation.py`):

The translation model uses **tRNA charging dynamics**:

```python
def evolve(self):
    # 1. Get ribosome positions and sequences being translated
    ribosome_data = self.active_ribosomes.molecules()
    
    # 2. Calculate amino acid demand from sequences
    sequence_elongations = self.proteinSequences[ribosome_data['protein_index']]
    
    # 3. tRNA charging kinetics (steady-state approximation)
    #    Synthetase + aa + tRNA_uncharged ⇌ tRNA_charged
    fraction_charged = calculate_trna_charging(...)
    
    # 4. Elongation rate depends on charged tRNA availability
    actual_rate = base_rate × fraction_charged
    
    # 5. Stochastic polymerization
    sequences_elongated, aas_consumed = polymerize(
        ribosome_data, actual_rate, aa_counts)
    
    # 6. Complete proteins: release ribosome, create protein
```

**Key biology captured:**
- Ribosome stalling when amino acids scarce
- Codon-dependent translation rates
- Competition between mRNAs for ribosomes

### Metabolism (Flux Balance Analysis)

**Metabolism** (`models/ecoli/processes/metabolism.py`) uses FBA:

```python
class FluxBalanceAnalysisModel:
    def solve(self, exchange_constraints, enzyme_counts):
        # 1. Build stoichiometric matrix S (metabolites × reactions)
        #    S·v = 0  (steady-state assumption)
        
        # 2. Set flux bounds from:
        #    - Enzyme counts (upper bound = kcat × [enzyme])
        #    - Nutrient uptake limits (from environment)
        #    - Reversibility constraints
        
        # 3. Objective: maximize biomass production
        #    max: c·v  where c is biomass reaction coefficients
        
        # 4. Solve linear program
        fluxes = scipy.optimize.linprog(...)
        
        # 5. Return metabolite production/consumption rates
```

**Integration with other processes:**
- Provides amino acids for translation
- Provides NTPs for transcription
- Consumes nutrients from environment
- Produces energy (ATP) for all processes

### DNA Replication (ChromosomeReplication)

```
Inputs:
  - Replication fork positions
  - DnaA-ATP concentration (initiation trigger)
  - dNTP availability

Algorithm:
  1. Check initiation conditions (DnaA-ATP > threshold)
  2. If initiating: create new replication forks at oriC
  3. Advance each fork by elongation_rate × timestep
  4. Track each fork's position independently
  5. When forks meet at terminus: chromosome complete

State tracked:
  - Number of chromosomes (1 → 2 during replication)
  - Fork positions (bp coordinates)
  - Replication progress (affects gene copy number)
```

### Cell Division

Division occurs when:
1. Chromosome replication complete (2 full chromosomes)
2. D-period elapsed (time for chromosome segregation)
3. Cell mass approximately doubled

```python
def evolveState(self):
    if (self.chromosome_complete and 
        self.d_period_complete and
        self.mass_doubled):
        self.trigger_division()
```

Division produces two daughter cells with:
- Half the cytoplasmic molecules (stochastic partitioning)
- One complete chromosome each
- Inherited state saved for next generation

---

## Outputs: What Data Comes Out

### Output Directory Structure

```
out/{sim_dir}/
├── kb/
│   └── simData.cPickle              # Fitted parameters (180 MB)
├── wildtype_000000/                 # Variant and index
│   └── 000000/                      # Random seed
│       └── generation_000000/
│           └── 000000/              # Daughter cell index
│               ├── simOut/          # Simulation data
│               │   ├── Main/        # Time, timestep
│               │   ├── Mass/        # Cell mass, growth
│               │   ├── BulkMolecules/  # All molecule counts
│               │   ├── UniqueMolecules/ # Individual molecules
│               │   ├── FBAResults/  # Metabolic fluxes
│               │   └── ...          # 18 listener directories
│               └── plotOut/         # Analysis figures
```

### Listener Data (18 Listeners)

Each listener records specific data each timestep:

| Listener | Data Recorded |
|----------|--------------|
| **Mass** | Total mass, dry mass, protein/RNA/DNA mass fractions |
| **BulkMolecules** | Counts of all ~3,000 bulk molecule species |
| **UniqueMolecules** | Individual RNAP, ribosome, replication fork states |
| **FBAResults** | Reaction fluxes, exchange rates, shadow prices |
| **GrowthLimits** | Which resources limit growth (AAs, NTPs, energy) |
| **RnaSynthProb** | Gene expression probabilities and actual synthesis |
| **RibosomeData** | Active ribosomes, stalling events, translation rates |
| **ReplicationData** | Fork positions, initiation events, chromosome count |
| **EnzymeKinetics** | Enzyme activities, metabolite concentrations |

### Reading Output Data

```python
from wholecell.io.tablereader import TableReader
import os

sim_out = "out/manual/wildtype_000000/000000/generation_000000/000000/simOut"

# Read time series
main = TableReader(os.path.join(sim_out, "Main"))
time = main.readColumn("time")  # seconds

# Read mass data
mass = TableReader(os.path.join(sim_out, "Mass"))
dry_mass = mass.readColumn("dryMass")  # femtograms
growth_rate = mass.readColumn("growth")  # fg/s

# Read molecule counts
bulk = TableReader(os.path.join(sim_out, "BulkMolecules"))
counts = bulk.readColumn("counts")  # shape: (timesteps, molecules)
molecule_ids = bulk.readAttribute("moleculeIds")  # list of IDs

# Find specific molecule
atp_idx = molecule_ids.index("ATP[c]")
atp_counts = counts[:, atp_idx]
```

---

## Visualization: Analysis Plots

### Analysis Organization

```
models/ecoli/analysis/
├── single/      # One cell, one generation (45+ plots)
├── multigen/    # One lineage across generations
├── cohort/      # Multiple seeds, statistical analysis
├── variant/     # Compare experimental conditions
└── parca/       # Parameter calculator validation
```

### Running Analysis

```bash
# Single-cell analysis (most detailed)
python runscripts/manual/analysisSingle.py out/manual

# Run specific plot
python runscripts/manual/analysisSingle.py out/manual --plot aaCounts

# Run plot group (defined by TAGS in __init__.py)
python runscripts/manual/analysisSingle.py out/manual --plot CORE
```

### Example Plots and What They Show

**aaCounts.py** - Amino acid dynamics:
- 21 subplots, one per amino acid
- Shows how AA pools fluctuate during growth
- Reveals metabolic bottlenecks

**growthDynamics.py** - Mass accumulation:
- Dry mass vs time (exponential growth)
- Protein/RNA/DNA mass fractions
- Comparison to expected doubling time

**rnaSynthProbs.py** - Gene expression:
- Probability of each gene being transcribed
- Comparison to expected steady-state levels
- Identifies dysregulated genes

**metabolicFluxes.py** - Metabolic state:
- Central carbon metabolism fluxes
- Amino acid biosynthesis rates
- Energy production (ATP flux)

---

## What Can You Use It For?

### Predicting Gene Knockout Phenotypes

```bash
# Run knockout variant
python runscripts/manual/runSim.py -v gene_knockout 0 2
```

The model predicts:
- Whether cell can grow (essential gene?)
- Growth rate changes
- Metabolic rewiring
- Expression compensation

### Testing Growth Conditions

```bash
# Different carbon source (variant index maps to media)
python runscripts/manual/runSim.py -v condition 1 2
```

Simulates growth in:
- Different carbon sources (glucose, acetate, glycerol)
- Amino acid supplementation
- Stress conditions (osmotic, oxidative)

### Understanding Regulation

The model captures:
- Transcription factor dynamics
- ppGpp stress response
- Two-component signaling systems
- Feedback loops in metabolism

### Generating Hypotheses

By comparing simulation to experiment:
- Identify missing regulation
- Find incorrectly parameterized processes
- Discover emergent behaviors

---

## Limitations and Caveats

### What's Not Modeled

- **Spatial organization** - Cell treated as well-mixed
- **Membrane processes** - No detailed transport kinetics
- **Cell shape** - No morphological changes
- **Stochastic gene expression** - Uses mean-field approximation for most regulation
- **Post-translational modifications** - Limited to essential modifications

### Parameter Uncertainty

- Many kinetic parameters estimated or fitted
- Some protein half-lives from different conditions
- Metabolic concentrations from mixed literature sources

### Computational Requirements

- **ParCa**: ~18 minutes, <8 GB RAM
- **Simulation**: ~10 minutes per cell cycle, <4 GB RAM
- **Full workflow**: 30+ minutes for one complete run

---

## Where to Look Next, by Goal

**To understand the simulation architecture:**
- `wholecell/sim/simulation.py`
- `models/ecoli/sim/simulation.py`
- `wholecell/processes/process.py`

**To change biology:**
- a relevant file in `models/ecoli/processes/`
- matching listeners in `models/ecoli/listeners/`
- any related reconstruction code in `reconstruction/ecoli/`

**To change fitted assumptions or conditions:**
- `reconstruction/ecoli/fit_sim_data_1.py`
- other files under `reconstruction/ecoli/`
- variants under `models/ecoli/sim/variants/`

**To add or debug plots:**
- `models/ecoli/analysis/analysisPlot.py`
- the relevant analysis subdirectory
- `wholecell/io/tablereader.py`

**To improve user workflows:**
- `runscripts/manual/`
- `wholecell/webapp/` (the Dash app — `app.py`, `jobs.py`, `results.py`,
  and the `tabs/` subdirectory)

---

## Quick Reference

### Running the Model

```bash
# 1. Compile Cython (required once)
make compile

# 2. Run parameter calculator (~18 min)
python runscripts/manual/runParca.py out/test

# 3. Run simulation (~10 min)
python runscripts/manual/runSim.py out/test

# 4. Generate plots
python runscripts/manual/analysisSingle.py out/test
```

### Key Files

| Purpose | File |
|---------|------|
| Parameter fitting | `reconstruction/ecoli/fit_sim_data_1.py` |
| Simulation loop | `wholecell/sim/simulation.py` |
| E. coli processes | `models/ecoli/processes/*.py` |
| E. coli listeners | `models/ecoli/listeners/*.py` |
| Analysis plots | `models/ecoli/analysis/single/*.py` |
| Raw data | `reconstruction/ecoli/flat/*.tsv` |

### Process Methods

Every process implements:
```python
def initialize(self, sim, sim_data):  # Setup
def calculateRequest(self):            # Declare needs
def evolveState(self):                 # Execute biology
```

### Molecule Types

| Type | Examples | State Container |
|------|----------|-----------------|
| Bulk | ATP, amino acids, proteins | `BulkMolecules` |
| Unique | Active RNAPs, ribosomes, replication forks | `UniqueMolecules` |
| External | Glucose, oxygen | `LocalEnvironment` |
