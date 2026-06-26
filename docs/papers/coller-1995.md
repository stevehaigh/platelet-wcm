# Coller (1995) — Blockade of platelet GPIIb/IIIa receptors as an antithrombotic strategy

- **Full citation:** Coller BS. *Blockade of Platelet GPIIb/IIIa Receptors as an Antithrombotic Strategy.* Circulation 92(9):2373–2380, 1995. DOI: 10.1161/01.CIR.92.9.2373.
- **Type:** review / perspective essay (bench-to-bedside narrative)
- **Local source:** source-info/calcium-papers/coller-1995-blockade-of-platelet-gpiib-iiia-receptors-as-an-antithrombotic-strategy.epub
- **Used in the model for:** biological/clinical grounding of the integrin αIIbβ3 (GPIIb/IIIa) activation module and its knockout (Glanzmann thrombasthenia / αIIbβ3 antagonist).

## What the paper does
Coller recounts the development of the monoclonal antibody 7E3 — from a laboratory reagent into the antithrombotic drug abciximab (c7E3 Fab / ReoPro, FDA-approved Dec 1994) — and lays out the rationale for choosing the platelet GPIIb/IIIa (integrin αIIbβ3) receptor as a therapeutic target, drawing on Glanzmann thrombasthenia, platelet–fibrinogen biology, the integrin family, and the RGD binding mechanism.

## Key findings / values the model uses
- **GPIIb/IIIa = integrin αIIbβ3**: GPIIb is the αIIb subunit, GPIIIa the β3 subunit, held in a **calcium-dependent complex**. It is platelet-specific and is the **final common pathway** for aggregation, regardless of agonist.
- **Copy number ~45 000 receptors/platelet** (Table 1; antibody-binding estimate ~40 000–80 000) — "probably the most dense adhesion/aggregation receptor on any cell."
- **Inside-out activation**: platelet activation drives a **conformational change** that switches GPIIb/IIIa to high affinity for fibrinogen and von Willebrand factor; multivalent ligands then crosslink platelets → aggregation.
- **Glanzmann thrombasthenia** = inherited GPIIb/IIIa deficiency → platelets fail to aggregate to *all* physiological agonists, despite intact adhesion; bleeding is mucocutaneous but spontaneous CNS haemorrhage is rare.
- **Dose-response of blockade**: <50% receptors available → significant inhibition; ~80% blockade → aggregation nearly abolished with only mild bleeding-time effect; 90% blockade → bleeding time 15–30 min.
- Aspirin blocks only arachidonic-acid metabolism → it is only a *partial* inhibitor; GPIIb/IIIa blockade is far more complete.

## Relevance to the platelet WCM
This is the clinical/biological reference behind `IntegrinActivation` (the αIIbβ3 inside-out resting⇌active switch; PAC-1 readout) and behind the two integrin knockout knobs: `integrin_act_scale=0` models an **αIIbβ3 antagonist / Glanzmann thrombasthenia** (no high-affinity integrin). The "final common pathway" framing motivates treating αIIbβ3 as the terminal output of the PKC/Rap1b cascade.

## Caveats / notes
This is a narrative review, not a kinetics source — it supplies the receptor identity, copy number, and the inside-out-activation concept, not rate constants. The model represents only the **per-cell affinity state** (PAC-1); **aggregation itself is inter-cellular and out of single-cell scope**. Fibrinogen-bound occupancy (a second possible readout) is not yet modelled. Text extracted from the EPUB and read directly.
