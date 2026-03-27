# What is wcEcoli and Why Am I Using It?

## What is wcEcoli?

wcEcoli is an open-source whole-cell computational model of *Escherichia coli*, developed by
the Covert Lab at Stanford University. It is the most complete whole-cell model of a
free-living organism currently available.

The project descends from the first-ever whole-cell model, published in Karr et al. (2012,
*Cell*), which simulated the entire life cycle of *Mycoplasma genitalium* — a minimal
bacterium with only ~470 genes. The *E. coli* model is substantially more complex: it covers
roughly 43% of the organism's characterised genes and integrates 17 distinct biological
processes running simultaneously.

**Key references:**

- Karr JR et al. (2012) *A whole-cell computational model predicts phenotype from genotype.*
  Cell 150(2):389–401. https://doi.org/10.1016/j.cell.2012.05.044
- Sun G, Ahn-Horst TA, Covert MW (2021) *The E. coli whole-cell modeling project.*
  EcoSal Plus 9(2). https://doi.org/10.1128/ecosalplus.ESP-0001-2020
- Ahn-Horst TA et al. (2022) *An expanded whole-cell model of E. coli links cellular physiology
  with mechanisms of growth rate control.* npj Syst Biol Appl 8:30.
  https://doi.org/10.1038/s41540-022-00242-9

---

## What does it simulate?

The model simulates a single *E. coli* cell from birth to division, advancing in one-second
timesteps. It tracks every major cellular process simultaneously:

- **Metabolism** — flux-balance analysis over the full metabolic network; glucose uptake,
  amino acid synthesis, energy production
- **Transcription** — RNA polymerase binding, transcript initiation and elongation, mRNA
  degradation
- **Translation** — ribosome assembly, polypeptide elongation, protein degradation
- **DNA replication** — replisome positioning, origin firing, chromosome segregation
- **Gene regulation** — transcription factor binding/unbinding, two-component signalling,
  ppGpp-mediated stringent response
- **Cell division** — mass-based division trigger, daughter cell initialisation

Each process uses the mathematical formalism best suited to it (stochastic for gene
regulation, ODE for tRNA charging, FBA for metabolism), and all share a common molecular
inventory that is partitioned and reconciled at every timestep.

Output is one-second-resolution time series for every tracked quantity: ~18 listener modules
recording mass fractions, metabolic fluxes, ribosome counts, replication fork positions,
amino acid levels, and more.

---

## Why is it scientifically significant?

Whole-cell models represent a qualitative shift from conventional pathway modelling. Rather
than asking "what does this pathway do in isolation?", a whole-cell model asks "what does the
cell do as a system?" This matters because:

- **Emergent behaviour.** The model reproduces the bacterial growth law (ribosome fraction
  scales linearly with growth rate) as an emergent consequence of its biochemistry, not as a
  built-in parameter. This kind of validation — predicting something that was not explicitly
  encoded — is the hallmark of a mechanistically complete model.
- **Genotype-to-phenotype.** Gene knockouts, nutrient shifts, and antibiotic perturbations
  can be simulated and their system-wide consequences observed, including unexpected
  off-target effects.
- **Integration across scales.** A single simulation spans molecular events (individual
  ribosome translation steps) and cell-level outcomes (doubling time, growth rate) without
  requiring separate models at each scale.

---

## Why am I using it?

I am using wcEcoli as the starting point for a proposed whole-cell model of the human
platelet. The rationale has two parts: scientific and practical.

### Scientific rationale

Platelets are clinically important but computationally understudied. Existing computational
models of platelet biology focus on individual pathways — calcium signalling, integrin
activation, coagulation — but no model integrates them into a single, self-consistent
simulation. A whole-cell model would allow questions such as:

- How does the kinetics of receptor activation propagate through calcium signalling to
  granule release?
- How does platelet energetics (ATP supply) constrain the speed or completeness of
  activation?
- What explains heterogeneity in activation thresholds across individual platelets?

These are precisely the kinds of cross-pathway, multi-scale questions that whole-cell
modelling is designed to address.

### Practical rationale

Building a whole-cell modelling framework from scratch is a multi-year engineering effort.
wcEcoli provides a mature, well-tested foundation:

- The codebase cleanly separates a **generic simulation engine** (`wholecell/`) from the
  **organism-specific model** (`models/ecoli/`). The engine — simulation loop, state
  management, I/O, listener framework — is organism-agnostic by design.
- The engine already handles the technical challenges that make whole-cell modelling hard:
  mixed mathematical formalisms within one simulation, discrete molecule tracking, parallel
  process execution with shared state, efficient binary I/O at one-second resolution.
- The web UI I have built works with any organism namespace, so a platelet model would
  immediately have a browser-based interface for running and inspecting simulations.

The practical cost of using wcEcoli is learning a large and complex codebase. The benefit is
avoiding years of framework engineering and starting with infrastructure that has already been
validated against decades of *E. coli* experimental data.

---

## Current status

I have the model running locally and on Azure. The web UI (on the `webapp` branch) is
functional and actively developed. I have conducted a detailed reusability analysis of every
wcEcoli component and documented the engine changes needed for a non-dividing cell type.

The next step is building a minimal platelet stub — a skeleton `models/platelet/` namespace
that runs through the existing engine, produces output visible in the web UI, and can be
incrementally extended with real biology.
