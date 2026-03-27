# The 4D Whole-Cell Model of JCVI-syn3A (MC4D)

**Source paper:** Thornburg ZR, Maytin A, Kwon J, Brier TA, Gilbert BR, Luthey-Schulten Z et al.
"Bringing the genetically minimal cell to life on a computer in 4D."
*Cell*, published 9 March 2026. DOI: [10.1016/j.cell.2026.02.009](https://doi.org/10.1016/j.cell.2026.02.009)

**Code repository:** [Luthey-Schulten-Lab/Minimal\_Cell\_4DWCM](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM)

---

## What it is

MC4D is a whole-cell computational model of **JCVI-syn3A**, a synthetically constructed
minimal bacterium with approximately 473 genes — one of the simplest organisms capable of
independent replication. The model simulates a complete ~100-minute cell cycle, capturing
every major biological process across three spatial dimensions and time (hence "4D").

This is the first whole-cell model to integrate full spatial resolution with a complete
biological description of cell growth and division. Where earlier whole-cell models
(including the original *M. genitalium* model and wcEcoli) treat the cell interior as
well-mixed compartments, MC4D tracks where molecules are and how they move.

---

## The organism: why JCVI-syn3A?

JCVI-syn3A is a stripped-down synthetic cell derived from *Mycoplasma mycoides* through
iterative gene removal, retaining only the minimal set of genes needed for growth and
division. With ~473 genes, it is dramatically simpler than *E. coli* (~4,300 genes). This
makes it tractable for a model that must account for every gene product — there are simply
fewer things to include.

The simplicity comes with a cost: roughly one third of syn3A's genes are of unknown
function, which the model must accommodate with placeholder or inferred behaviour.

---

## The four computational methods

The model is named "4D" both for its spatial dimensions and for the four distinct
mathematical formalisms it integrates simultaneously.

### 1. RDME — Reaction-Diffusion Master Equation
*Implemented via: Lattice Microbes (GPU-accelerated)*

The cell interior is discretised onto a 3D lattice. Molecules are represented as particles
that diffuse stochastically between lattice voxels and react with each other based on local
concentrations. This captures spatial heterogeneity: the fact that two molecules can only
react if they are in the same place at the same time.

Used for: the bulk of intracellular chemistry — protein synthesis (translation), molecular
diffusion, bimolecular reactions, ribosome dynamics with excluded-volume constraints.

### 2. CME — Chemical Master Equation
*Implemented as global stochastic reactions*

A well-mixed stochastic formalism for reactions that are best treated globally rather than
spatially. Used for tRNA charging reactions and transcription events, where the assumption
of spatial homogeneity is reasonable and the stochastic discrete-molecule treatment matters
(low copy numbers).

### 3. ODE — Ordinary Differential Equations
*Implemented via: odecell framework*

Deterministic kinetic integration of the metabolic network. Enzyme concentrations are
extracted from the RDME simulation and passed into the ODE solver, so metabolic fluxes
reflect the actual current protein composition of the cell. Outputs include dNTP pool
concentrations, which are fed back to control DNA replication rates.

### 4. BD — Brownian Dynamics
*Implemented via: btree\_chromo with Kokkos-enabled LAMMPS*

Explicit physical simulation of chromosomal DNA as a polymer. The chromosome is represented
as a chain of beads subject to thermal fluctuations, constrained to the RDME lattice.
Structural maintenance of chromosome (SMC) proteins (analogous to condensin) and
topoisomerases control chromosome compaction and segregation. This is the most
computationally demanding component and requires the LAMMPS molecular dynamics engine.

---

## How the four methods communicate

Integration is managed by a central Python module, `Hook.py`. This module periodically
interrupts the RDME simulation to:

1. Extract current enzyme and molecular counts from the RDME lattice
2. Pass enzyme concentrations to the ODE solver; run metabolic integration
3. Pass transcription events and tRNA status to the CME solver; advance gene expression
4. Synchronise chromosome state with the BD simulation
5. Write updated molecular counts back to the RDME lattice

This coupling means that metabolism, gene expression, and chromosome dynamics all influence
each other dynamically throughout the simulation, rather than being run in isolation.

---

## Biological processes covered

| Process | Method | Key detail |
|---|---|---|
| Metabolism | ODE | Kinetic network; dNTP output constrains replication |
| Transcription | CME | Global stochastic; mRNA production and degradation |
| Translation | RDME | Spatially resolved; ribosome excluded-volume |
| tRNA charging | CME | Global stochastic |
| DNA replication | BD + ODE | Chromosome polymer dynamics; rate set by dNTP pools |
| Chromosome segregation | BD | SMC and topoisomerase proteins modelled explicitly |
| Membrane growth | RDME | Lipid and membrane protein insertion drives morphology |
| Cell division | RDME + BD | Triggered by chromosome segregation completion |
| Molecular diffusion | RDME | All major cytoplasmic species |

---

## Spatial resolution and scale

- **Lattice**: 3D voxel grid representing the cell interior
- **Cell morphology**: initialised as a sphere; morphological change during growth constrained
  by fluorescence imaging data from experiments
- **Timescale**: full ~100-minute cell cycle
- **Stochasticity**: because of the stochastic formalisms (RDME, CME), each simulation
  replicate produces a unique cell — population heterogeneity emerges naturally

---

## Validation

The model is validated against a wide array of experimental measurements:

- **Doubling time** — reproduced quantitatively
- **Origin-to-terminus ratio** — validated against DNA sequencing data
- **mRNA half-lives** — matched to experimental measurements
- **Protein distributions** — spatial and abundance patterns compared to fluorescence data
- **Ribosome counts** — quantitative agreement
- **Daughter cell heterogeneity** — stochastic partitioning of molecules matches
  experimentally observed variability between siblings

The authors emphasise that "assimilation of a wide array of experiments is necessary for
construction and validation" — the model is tightly constrained by data at multiple levels.

---

## Software and technical requirements

| Component | Software | Notes |
|---|---|---|
| RDME solver | Lattice Microbes | GPU required |
| ODE integrator | odecell | Python-based |
| Chromosome dynamics | btree\_chromo + LAMMPS | Kokkos-enabled build required |
| Chromosome initialisation | sc\_chain\_generation | Generates initial polymer configs |
| Membrane (optional) | FreeDTS | For membrane shape analysis |
| Orchestration | Hook.py (custom) | Central integration module |

Running a full simulation requires a GPU (for Lattice Microbes) and a correctly configured
LAMMPS build with Kokkos support. The model runs as a single Python entry point
(`Whole_Cell_Minimal_Cell.py`) with parameters for simulation duration, output directory,
GPU device selection, and random seeds for reproducibility.

Input data lives in `input_data/` and includes kinetic parameters (xlsx), genomic sequences
(GenBank), metabolic network definitions (SBML), and DNA polymer properties.

---

## Key scientific contribution

The central claim of the paper is that spatial heterogeneity within the cell — the fact that
molecules are not uniformly distributed — materially affects the biochemical reactions that
control cellular phenotype. By resolving this spatial structure, MC4D can predict behaviours
that well-mixed models cannot, including the spatial patterns of protein synthesis, the
physical dynamics of chromosome segregation, and the sources of cell-to-cell variability
in division outcomes.

This is a significant step beyond the original *M. genitalium* whole-cell model (which used
the same organism's precursor but had no spatial resolution) and beyond wcEcoli (which
achieves greater biological coverage but treats the cell as a set of well-mixed
compartments).

---

## Limitations and context

- Requires specialist software dependencies and GPU hardware — harder to run than wcEcoli
- The organism (JCVI-syn3A) is simpler than any naturally occurring bacterium; many results
  may not generalise directly
- ~1/3 of syn3A genes have unknown function; the model must make assumptions about these
- No public web interface or high-level user tooling; the model is a research prototype
- Computational cost is substantially higher than non-spatial whole-cell models
