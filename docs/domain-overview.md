# Domain overview — what this project is and why

> Orientation doc for AI assistants and humans. For the *code* layout see
> [`codebase-overview.md`](codebase-overview.md); for the *engine* design see
> [`architecture.md`](architecture.md); for the source literature see
> [`papers/`](papers/).

## The one-paragraph version

This is a **whole-cell model (WCM) of a human platelet**, built as an MRes
dissertation project. It simulates the biochemistry a platelet runs when it
switches from quiescent to activated: agonists bind surface receptors, a Gq
cascade produces IP₃, IP₃ releases Ca²⁺ from the internal store, and the Ca²⁺
signal drives the functional outputs that make a platelet useful — granule
secretion, thromboxane production, and integrin (αIIbβ3) activation. The
simulation engine is reused from [CovertLab/wcEcoli](https://github.com/CovertLab/wcEcoli)
(the *E. coli* whole-cell model); **all of the E. coli biology was removed** and
replaced with platelet-specific processes, listeners, and a hand-built parameter
set. Current scope is a **single resting→activated platelet over seconds–minutes**
(no division, no aggregation).

## The biological problem

Platelets are the small, anucleate cell fragments that stop bleeding. When a
vessel is injured they:

1. **Sense** the injury via surface receptors (collagen via GPVI; thrombin via
   PAR1/PAR4; ADP via P2Y1/P2Y12; ATP via P2X1; thromboxane via TP).
2. **Signal** internally — the central integrator is a rise in cytosolic
   **Ca²⁺**, produced by the GPCR → Gq → PLCβ → IP₃ → IP₃R cascade releasing the
   internal Ca²⁺ store (the **Dense Tubular System**, the platelet's ER), backed
   up by **store-operated Ca²⁺ entry (SOCE)** from outside the cell.
3. **Respond** — Ca²⁺ (with PKC) drives shape change, **dense/α-granule
   secretion** (which releases more ADP → autocrine amplification),
   **thromboxane A₂ synthesis** (another autocrine amplifier), and **inside-out
   activation of integrin αIIbβ3** (the receptor that binds fibrinogen and lets
   platelets aggregate).
4. **Are restrained** by an inhibitory axis — endothelial prostacyclin (PGI₂) and
   nitric oxide raise cAMP/cGMP, activating PKA/PKG, which brake the activation
   machinery. Antiplatelet drugs (aspirin, clopidogrel, etc.) act on these
   pathways.

Why platelets are an interesting WCM target: they are **anucleate**, so there is
no transcription/translation/cell-cycle machinery to model — the system is a
**closed signalling + secretion network** running off a pre-loaded protein and
metabolite inventory. That makes a "whole-cell" scope tractable in a way it is
not for a dividing, transcribing cell.

## How we tackle it

- **Reuse the engine, not the ontology.** wcEcoli gives a tested
  process/state/listener architecture, discrete time-stepping, resource
  partitioning, columnar I/O, and analysis plumbing. We keep that and write new
  *platelet* processes inside it. See [`architecture.md`](architecture.md).
- **Build biology incrementally, calcium-first.** The Ca²⁺ core was built and
  validated first (against Dolan & Diamond 2014), then each downstream module
  (PKC feedback, secretion, thromboxane, integrin) and the inhibitory axis
  (P2Y12/Gi/cAMP/PKA) were layered on, each normalised so the resting state and
  the validated Ca²⁺ transient are preserved.
- **Parameters are data, not code.** Rate constants and the molecule inventory
  live in versioned TOML/TSV files under `reports/params/` with per-value
  citations, loaded at import. You change biology by editing data.
- **Ground every number in a paper.** Numeric values are traced to a primary
  source and verified against the PDF (`reports/data/calcium-data-provenance.md`,
  and the per-paper notes in [`papers/`](papers/)). Findings are checked for
  *direction*, not just value.

## What the model does and does not cover

**In scope (built):**

- Ca²⁺ homeostasis & mobilisation: IP₃R (Li-Rinzel/Sneyd-Dufour), SERCA (Dode),
  PMCA (Caride CaM-coupled), SOCE (Hoover-Lewis/Dolan MWC), CaM buffering,
  multi-buffer DTS store, MCU mitochondrial uptake.
- GPCR cascade: P2Y1, PAR1, PAR4, P2X1, TP, P2Y12 → Gq/Gi → PLCβ → PI cycle →
  endogenous IP₃ (not a forced curve).
- PKC feedback (P2Y1 desensitisation + PLCβ phosphorylation).
- Functional outputs: granule secretion + autocrine ADP, thromboxane synthesis +
  autocrine TP loop, integrin αIIbβ3 inside-out activation (PAC-1 readout) via a
  PI3K/Akt→Rap1b arm.
- Inhibitory axis: P2Y12/Gi/cAMP/PKA/VASP, with cAMP-raising agents (PGI₂,
  forskolin, cilostazol).
- Pharmacology knobs: aspirin (COX-1), clopidogrel (P2Y12), MCU inhibitor,
  αIIbβ3 antagonist / Glanzmann KO, and assorted scale/knockout knobs.

**Out of scope (not built / future):**

- Cell division, transcription, translation (anucleate — deliberately absent).
- Aggregation and any multi-cell behaviour (αIIbβ3 *occupancy* by fibrinogen is
  not yet modelled; only the per-cell affinity state / PAC-1 readout is).
- GPVI → Syk → PLCγ2 collagen arm (a candidate future pathway; see Dunster 2015).
- The kinetics-as-data scaffold is **calcium-only**; other pathways carry their
  constants in Python dataclasses.
- Spatial / microdomain detail; stochastic single-cell variability (the model is
  deterministic ODEs — see the "stochastic bottleneck" note in the provenance
  doc).

## The headline validation, and its honest limits

The model's original validation target is **Dolan & Diamond 2014, Fig. 4** — the
cytosolic Ca²⁺ transient with and without extracellular Ca²⁺ ("Dolan 5/5"). As
the model grew past being *just* the Ca²⁺ transient, that single figure stopped
being a sufficient measure of correctness. We now treat 5/5 as a **regression
invariant** on the Ca²⁺ core and judge the rest of the model on a **portfolio**
of subsystem-specific assays. This is discussed candidly — including why the Ca²⁺
signal is "clamped" and what that means for what the model can predict — in
[`validation-and-regressions.md`](validation-and-regressions.md).
