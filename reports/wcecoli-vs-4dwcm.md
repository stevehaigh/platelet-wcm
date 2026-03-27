# wcEcoli vs MC4D: A Comparison of Two Whole-Cell Modelling Approaches

These are the two most complete whole-cell models currently available. They share a common
ambition — simulate an entire cell — but make substantially different choices about organism,
formalism, and what "complete" means. Understanding the differences matters for deciding
which approach, or which combination, is most useful for extending to human platelets.

---

## Side-by-side overview

| Dimension | wcEcoli | MC4D |
|---|---|---|
| **Organism** | *E. coli* (wild-type, ~4,300 genes) | JCVI-syn3A (~473 genes, synthetic minimal) |
| **Gene coverage** | ~43% of characterised genes | ~100% (all ~473 genes, including unknowns) |
| **Spatial resolution** | None — well-mixed compartments | Full 3D (RDME lattice + BD chromosome) |
| **Cell cycle** | ~25–130 min depending on condition | ~100 min |
| **Metabolism** | Flux-balance analysis (FBA) | ODE kinetic network |
| **Transcription** | Stochastic (CME-like), well-mixed | CME, global stochastic |
| **Translation** | Stochastic, well-mixed | RDME, spatially resolved |
| **Chromosome** | Discrete replisome tracking (position only) | Brownian dynamics polymer (LAMMPS) |
| **Formalisms** | FBA + stochastic + ODE | RDME + CME + ODE + BD |
| **Stochasticity** | Partial (gene regulation, translation) | Full (all molecular processes) |
| **GPU required** | No | Yes (RDME via Lattice Microbes) |
| **Web interface** | Yes (custom Dash app) | None |
| **Framework reuse** | Explicitly designed for it | Research prototype |
| **Published** | npj Sys Bio Appl, 2022 | Cell, March 2026 |
| **Lab** | Covert Lab, Stanford | Luthey-Schulten Lab, UIUC |

---

## Organism choice: complexity vs completeness

This is the most fundamental difference. wcEcoli models *E. coli* — a well-characterised,
non-minimal organism with thousands of genes and decades of experimental data. The trade-off
is that the model can only cover ~43% of genes; the rest are either not yet modelled or not
characterised well enough to include.

MC4D models JCVI-syn3A, a synthetic cell engineered to be as simple as possible. With ~473
genes, it is feasible — in principle — to include every gene product. The trade-off is that
roughly a third of those genes have unknown function, and the organism is too stripped-down
to exhibit many of the regulatory behaviours (stress responses, growth-rate adaptation) that
make *E. coli* biologically interesting.

**For platelet modelling:** Neither organism is a direct analogue, but the wcEcoli philosophy
— deep coverage of a complex, well-studied cell — is more relevant. Platelets are well
characterised biochemically; the goal is integrating known biology, not filling in unknowns.

---

## Spatial resolution: the central architectural divide

wcEcoli treats the cell as a set of well-mixed compartments (cytoplasm, periplasm,
extracellular). Molecules exist as integer counts, not positions. This is computationally
efficient and sufficient for many questions, but it cannot address anything that depends on
*where* a molecule is.

MC4D places every molecule on a 3D lattice and simulates diffusion explicitly. This reveals
spatial heterogeneity — gradients of signalling molecules, localised translation near the
membrane, chromosome-excluded cytoplasmic volumes — that well-mixed models cannot see.
The paper argues these spatial effects materially influence phenotype.

The cost is substantial: spatial simulation at full resolution requires GPU acceleration and
is computationally much more demanding. A full cell cycle in MC4D requires specialist
hardware; wcEcoli runs on a laptop in ~10 minutes.

**For platelet modelling:** Spatial resolution could matter significantly. Platelets have
highly organised internal compartments (dense granules, alpha-granules, open canalicular
system, mitochondria), and calcium signalling involves spatial waves propagating from
stores to the plasma membrane. A spatially resolved approach could model these dynamics
more faithfully. However, the added complexity and hardware requirements are a significant
practical barrier for initial modelling work.

---

## Metabolic treatment: FBA vs ODE kinetics

wcEcoli uses flux-balance analysis for metabolism: it finds the optimal flux distribution
through the metabolic network at each timestep, given the current enzyme activities and
nutrient availability. FBA does not require kinetic rate constants for every reaction, which
is an advantage when those constants are unknown. The downside is that FBA assumes
steady-state — it cannot model transient metabolic dynamics.

MC4D uses ODE kinetic modelling for metabolism: reactions have explicit rate laws, and
metabolite concentrations evolve over time. This captures transient dynamics (e.g., rapid
dNTP depletion during replication) but requires kinetic parameters for every reaction.

**For platelet modelling:** Transient dynamics matter enormously for platelets. Activation
happens over seconds to minutes, and the ATP/ADP ratio, calcium concentration, and signalling
second messengers all change rapidly. ODE kinetics is the more appropriate formalism for
platelet metabolism, and MC4D's approach is better aligned with the biology.

---

## Chromosome and replication treatment

In wcEcoli, DNA replication is modelled as discrete replisomes advancing along the
chromosome at rates determined by nucleotide availability and polymerase speed. Position is
tracked as a single coordinate. There is no physical model of the chromosome as a polymer.

MC4D models the chromosome as a Brownian dynamics polymer using LAMMPS molecular dynamics.
SMC proteins (condensin analogues) and topoisomerases are represented explicitly and control
compaction and segregation. dNTP pools from the metabolic ODE regulate replication speed,
creating feedback between metabolism and chromosome dynamics.

**For platelet modelling:** This difference is largely irrelevant — platelets have no
nucleus and no DNA replication. Neither model's chromosome treatment will be used.

---

## Stochasticity

Both models are stochastic, but MC4D is more thoroughly so. In wcEcoli, some processes
(gene regulation, molecule partitioning) are stochastic; others (metabolism via FBA,
initiation rules) are deterministic. In MC4D, all molecular processes are governed by
stochastic formalisms (RDME or CME), and the result is that every simulated replicate
produces a genuinely unique cell — population heterogeneity is an intrinsic output of
the model, not a post-hoc analysis.

This matters for platelet biology. Platelets show substantial cell-to-cell heterogeneity in
activation thresholds, granule content, and response kinetics. A fully stochastic model
would capture this naturally.

---

## Framework design and reusability

wcEcoli is explicitly architected for reuse. The separation between `wholecell/` (generic
engine) and `models/ecoli/` (organism-specific biology) is a deliberate design choice. The
framework has already been extended in the literature (colony-level simulations in Skalnik
et al. 2023) and is designed so a new organism namespace can be added.

MC4D is a research prototype. The code is open-source and well-documented, but it is not
architected as a reusable framework. Using it as a starting point for a different organism
would require substantial reworking of the core infrastructure.

---

## Tooling and accessibility

wcEcoli ships with 80+ analysis scripts, a custom binary data format with reader/writer
libraries, and (via this project) a browser-based web UI for running and inspecting
simulations without command-line expertise.

MC4D has no equivalent tooling. It produces trajectory files that require custom analysis.
There is no web interface, no preset simulation runner, and no built-in analysis scripts.
Running it requires GPU hardware, a custom LAMMPS build, and a conda environment with five
specialist dependencies.

---

## What each model is best suited for

| Question | Better model |
|---|---|
| How does growth rate depend on nutrient availability? | wcEcoli |
| How does gene expression vary across the cell spatially? | MC4D |
| What is the systems-level effect of a gene knockout? | wcEcoli |
| How does chromosome segregation affect daughter-cell heterogeneity? | MC4D |
| How do signalling cascades interact with metabolism? | wcEcoli (closer) |
| What is the spatial source of cell-to-cell variability? | MC4D |
| Running many simulations quickly for statistical analysis | wcEcoli |
| First principles spatial modelling of organelle dynamics | MC4D |

---

## Implications for the platelet model

A platelet whole-cell model sits somewhere between the two approaches in its requirements.
Some considerations:

- **Spatial heterogeneity matters.** Calcium waves, granule localisation, and membrane
  receptor clustering are inherently spatial. MC4D's RDME approach would handle these more
  faithfully than wcEcoli's well-mixed compartments.
- **Transient dynamics matter.** Platelet activation happens over seconds; FBA's
  steady-state assumption is unsuitable. ODE kinetics (as in MC4D) is the right formalism
  for the metabolic and signalling processes.
- **No chromosome needed.** The MC4D chromosome machinery (BD + LAMMPS) is not relevant to
  a platelet model and would not be carried over.
- **Framework reusability is practical.** Starting from wcEcoli's generic engine is more
  practical than reworking MC4D for a different organism. The key question is whether
  spatial resolution can be added later, once a working non-spatial stub exists.
- **A pragmatic path:** build the initial platelet model using the wcEcoli framework (fast,
  reusable, well-tooled), using ODE kinetics for calcium and signalling rather than FBA.
  If spatial effects prove important, the RDME approach from MC4D could be adopted for
  a second-generation model — or the two frameworks could be consulted for methodological
  guidance even if not directly reused.
