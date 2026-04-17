# Zotero Literature Analysis for Platelet Whole-Cell Model

*Generated: 2026-04-12. Covers 4 Zotero collections: Platelet Biology (55), Platelet Modelling (67), platelet_new_refs (16), Cell Modelling (18). Many items overlap across collections.*

---

## Summary

Of ~90 unique papers across the four collections, ~25 are already identified as key sources. This analysis focuses on **NEW papers not previously catalogued** that offer actionable contributions: quantitative parameters, implementable ODE models, or validation data.

---

## 1. Calcium Dynamics (Phase 3a priority)

### HIGH PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Kleppe 2018** — "Mathematical Modelling of NO/cGMP/cAMP Signalling in Platelets" | **ODE model of the inhibitory NO/cGMP/cAMP pathway** including PDE2, PDE3, PDE5, PKA, PKG. This is the main *brake* on platelet activation and is absent from all your current Ca2+ sources. Rate constants provided. Essential for modelling sub-threshold responses and aspirin-like inhibition. |
| **Fernández 2021** — "Platelet calcium signaling by GPCR and ITAM-linked receptors regulating anoctamin-6 and procoagulant activity" | Quantitative Ca2+ traces for GPCR vs ITAM pathways; defines the Ca2+ threshold for TMEM16F/anoctamin-6 activation (PS exposure). Provides the **decision criterion** for procoagulant vs aggregatory platelet commitment. |
| **Ghatge 2026** — "The mitochondrial calcium uniporter regulates calcium dynamics to drive platelet function, bioenergetics, and thrombosis" | Most recent paper on **MCU (mitochondrial calcium uniporter)** — quantifies how mito Ca2+ uptake shapes cytoplasmic Ca2+ transients and couples to aerobic glycolysis. Directly extends Dolan & Diamond's 34-species model with a mitochondrial compartment. |
| **Shehwar 2025** — "Platelets and mitochondria: the calcium connection" | Review synthesising mitochondrial Ca2+ uptake/release kinetics, MPTP opening thresholds, and the link to procoagulant subpopulation formation. Provides parameter ranges for MCU, NCLX, MPTP. |
| **Balabin (undated)** — "Personalization of a computational systems biology model of blood platelet calcium signaling" | Extends Sveshnikova group's Ca2+ model with **patient-specific parameterisation**. Shows which parameters vary most between individuals — useful for sensitivity analysis and identifying which rate constants matter most. |

### MEDIUM PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Ajanel 2025** — "Mitochondrial Calcium Uniporter Regulates ITAM-Dependent Platelet Activation" | Experimental data on MCU role specifically in GPVI/ITAM pathway (complements Ghatge 2026 which focuses on thrombin/GPCR). |
| **Tantiwong et al. 2025** — "Extended Modelling of Molecular Calcium Signalling in Platelets by RNN and PLS" | Neural network fits to Ca2+ curves under combinatorial agonist + inhibitor conditions. Not directly implementable as ODEs, but the **fitted Ca2+ trace dataset** across multiple conditions is excellent validation data. |
| **Shankar et al. 2025** — "Multiscale simulations with patient-specific NN models of platelet Ca2+ signaling" | Neural-network Ca2+ model trained on 10 donors × 6 agonists. Shows inter-donor variability ranges. Useful as **validation benchmark** — your ODE model should reproduce similar Ca2+ mobilisation patterns. |

---

## 2. Receptor Signalling (Phase 3c)

### HIGH PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Dunster 2015** — "Regulation of Early Steps of GPVI Signal Transduction by Phosphatases: A Systems Biology Approach" | **Data-driven ODE model of GPVI → Syk → PLCγ2 pathway** including SHP-1/SHP-2 phosphatase regulation. Provides rate constants for the entire GPVI signalling cascade upstream of IP3. Fills the gap between receptor binding and the Purvis 2008 IP3 model. |
| **Mazet 2020** — "A model of the PI cycle reveals regulating roles of lipid-binding proteins" | **ODE model of the phosphatidylinositol (PI) cycle** — PIP2 → IP3 + DAG, including lipid-binding protein regulation. Provides rate constants for PI4K, PIP5K, PLC, and lipid transfer proteins. Bridges receptor activation to IP3 generation. |
| **Kleppe 2018** (also listed under Ca2+) | The cGMP/cAMP model includes **P2Y12 → Gi → AC inhibition** pathway — the main ADP inhibitory receptor mechanism. |

### MEDIUM PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Guidetti 2015** — "PI3K/Akt in platelet integrin signaling and implications in thrombosis" | Review of PI3K/Akt pathway in inside-out and outside-in integrin signalling. Qualitative pathway maps useful for model architecture, not quantitative parameters. |
| **Bye 2016** — "Platelet signaling: a complex interplay between inhibitory and activatory networks" | Comprehensive signalling network map showing cross-talk between pathways. Useful reference for ensuring model connectivity is complete. |
| **Hashemzadeh 2023** — "The ten main platelet receptors" | Covers PAR1, PAR4, P2Y1, P2Y12, TP, GPVI, GP Ib-IX-V, αIIbβ3, CLEC-2, PEAR1. Receptor copy numbers and affinities where available. |
| **Swieringa 2025** — "Platelet activation and signaling in thrombus formation" | Most recent comprehensive signalling review; useful for cross-checking pathway completeness. |

---

## 3. Granule Exocytosis (Phase 3b)

### HIGH PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Fitch-Tewfik & Flaumenhaft 2013** — "Platelet Granule Exocytosis: A Comparison with Chromaffin Cells" | Detailed SNARE machinery (VAMP-3/8, syntaxin-11, SNAP-23) and **dynamin-mediated fusion pore dynamics**. Extends Flaumenhaft 2011 with mechanistic detail on fusion pore opening/closing kinetics. |
| **Lopez et al. 2015** — "Relationship between calcium mobilization and platelet α- and δ-granule secretion" | Quantitative flow cytometry data showing **α-granule secretion precedes δ-granule secretion** with different Ca2+ dependencies. TRPC6 antibody blocks δ- but not α-granule release. Provides the Ca2+ threshold data needed for granule release functions. |
| **Kulkarni et al. 2021** — "Mitochondrial ATP generation essential for granule secretion but dispensable for aggregation" | Key constraint: **mito-ATP is required for secretion but not aggregation/procoagulant activity**. Directly informs how to couple the metabolism and granule processes. |

### MEDIUM PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **van Nispen tot Pannerden 2010** — "The platelet interior revisited: electron tomography" | 3D structural data on α-granule subtypes (tubular vs spherical) and membrane organisation. Useful for compartment volume estimates. |
| **Ambrosio 2025** — "The winding road to platelet α-granules" | α-granule biogenesis and NBEAL2 role. More relevant to megakaryocyte biology than platelet activation, but useful for understanding initial granule state. |
| **Gerda 2025** — "Experimental and Mathematical Model of Platelet Hemostasis Kinetics" | Mathematical model of platelet phenotypic transitions during hemostasis. May provide transition rate constants between resting → activated → procoagulant states. |

---

## 4. Metabolism (Phase 3d)

### HIGH PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Thomas et al. 2014** — "Network reconstruction of platelet metabolism" (Palsson group) | **Constraint-based metabolic reconstruction** of the human platelet (iAT-PLT-636, 636 reactions). Directly implementable as a metabolic process. Identifies aspirin resistance metabolic signatures. |
| **Bordbar et al. 2011** — "iAB-RBC-283: erythrocyte metabolism" | Proteomically-derived metabolic reconstruction of the **erythrocyte** (281 enzymes). Methodological template for anucleate cell metabolism. Many enzymes shared with platelets. |
| **Flora et al. 2023** — "Mitochondrial PDK2/4 contribute to platelet function by regulating aerobic glycolysis" | Quantifies the **PDK→PDH phosphorylation switch** that diverts pyruvate from OXPHOS to aerobic glycolysis upon activation. Provides glycoPER measurements and the molecular mechanism for metabolic switching (extends Aibibula 2018). |

### MEDIUM PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Paglia 2014** — "Comprehensive metabolomic study of platelets" | Metabolomics during storage — maps discrete metabolic phenotypes. Useful for validating resting-state metabolite concentrations. |
| **Melchinger 2019** — "Role of Platelet Mitochondria: Life in a Nucleus-Free Zone" | Review of mitochondrial function in anucleate context: OXPHOS, apoptosis signalling, bioenergetics. Conceptual framework for the metabolism process. |
| **Mazet 2018** — "Precise Quantification of Platelet Proteins and Their Phosphorylation States" | Methods paper for quantitative Western blotting of platelet proteins. Provides absolute copy numbers for key signalling proteins by immunoblotting (cross-validates Burkhart proteomics). |

---

## 5. Whole-Cell Modelling Methodology

### HIGH PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Thornburg 2022/2026** — JCVI-syn3A minimal cell WCM | Most advanced whole-cell kinetic model (493-gene minimal cell). The 2026 paper adds **4D spatial modelling**. Directly relevant methodology for integrating ODE subsystems within a discrete-timestep framework. The minimal cell is the closest analogue to a platelet (simple, few processes). |
| **Elsemman 2022** — "Yeast whole-cell model with compartment-specific proteome constraints" | **Proteome-constrained metabolic modelling** — shows how to use proteome data (like Burkhart) as constraints on metabolic flux. Directly applicable methodology. |

### MEDIUM PRIORITY

| Paper | What it contributes NEW |
|-------|------------------------|
| **Karr 2015** — "Principles of whole-cell modeling" | Methodological review by the wcEcoli architect. Framework design principles you're already using, but useful reference for dissertation methodology section. |
| **Szigeti 2017** — "Blueprint for human whole-cell modeling" | Roadmap for human cell WCMs. Discusses challenges specific to human cells that may apply to the platelet model. |
| **Georgouli 2023** — "Multi-scale models of whole cells" | Reviews progress and challenges in multi-scale WCMs. Good for dissertation literature review. |

---

## 6. Multiscale / Thrombus Context (Lower priority for MVP)

| Paper | Relevance |
|-------|-----------|
| **Shankar et al. 2022/2023/2025** (Diamond group) | 3D multiscale thrombus models. Your intracellular model could eventually plug into these frameworks. Not needed for MVP but important for positioning the work. |
| **Sorensen 1999** | Early platelet deposition/activation computational model. Historical reference. |
| **Chung et al. 2023** | QSP coagulation cascade review — identifies reusable model components. |
| **Diamond 2016** | Systems analysis overview of thrombus formation. |
| **Wang 2023** | Shear-mediated platelet adhesion dynamics. Relevant if modelling flow effects. |

---

## 7. Papers with Low Relevance to Current Model

The following categories can be deprioritised:
- **Platelet-cancer interaction** (~10 papers): Genitoni 2025, Anvari 2021, Ortiz-Otero 2020, Zuo, Hinterleitner 2021, Fabricius 2021, Ward, Eslami-S 2023, Tesakov 2025, Menter 2014
- **Megakaryocyte biology**: Italiano 2003/2017, McArthur 2018, Furniss 2024
- **Clinical/transfusion**: Kelly 2017, Gorog 2025
- **General reviews without quantitative data**: Chaudhary 2022, Franco 2015, Boilard 2025, Sharma 2022
- **Unrelated**: Silen 1975, Lefkowitz 1975, Tohei 1975, Sencanski 2014, Vascular diseases WHO page

---

## Recommended Actions

### Immediate (for Phase 3a — Calcium Dynamics)

- [ ] **Read Kleppe 2018 in full** — the NO/cGMP/cAMP ODE model fills a critical gap (inhibitory signalling)
- [ ] **Read Ghatge 2026 + Shehwar 2025** — add mitochondrial Ca2+ compartment to Dolan & Diamond model
- [ ] **Extract validation data from Fernández 2021** — Ca2+ traces for GPCR vs ITAM, TMEM16F threshold

### Next (for Phase 3b-c — Granules & Receptors)

- [ ] **Read Dunster 2015** — GPVI signalling ODE with rate constants
- [ ] **Read Mazet 2020** — PI cycle ODE connecting receptors to IP3
- [ ] **Read Lopez et al. 2015** — Ca2+ thresholds for α vs δ granule release

### For Metabolism Process

- [ ] **Evaluate Thomas 2014 iAT-PLT-636** — may be directly usable as the metabolism process
- [ ] **Read Flora 2023** — PDK/PDH switch mechanism for activation-state metabolism

### For Dissertation Methodology

- [ ] **Cite Thornburg 2022/2026** as closest methodological precedent for minimal/anucleate WCM
- [ ] **Cite Elsemman 2022** for proteome-constrained approach
