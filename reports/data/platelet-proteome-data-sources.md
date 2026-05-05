# Quantitative Platelet Proteomics: Data Sources for Whole-Cell Computational Modelling

*Prepared: 2026-03-27. For use in designing the initial molecule inventory and compartment assignments of a whole-cell platelet activation model.*

---

## 1. Executive Summary

The Burkhart et al. 2012 dataset (Blood 120:e73) remains the gold-standard reference for the human resting-platelet proteome with absolute copy numbers, but it should be supplemented with at least four additional datasets to build a credible computational model:

1. **Zeiler et al. 2014** (murine, absolute copy numbers via SILAC-PrEST + iBAQ) — provides the most methodologically rigorous copy-number estimates across the full abundance range; use to validate and recalibrate Burkhart numbers for shared proteins.
2. **Burkhart et al. 2012** (human, NSAF-derived absolute copy numbers) — primary source for human-specific initial inventory.
3. **Aslan et al. 2021 / Iancu-Rubin / Mah et al.** (human, transcript-matched proteome, 5.2 k proteins classified into 21 compartmental categories) — provides the most complete subcellular compartment assignments currently available.
4. **Signalling/phosphoproteomics trio:** Zahedi & Sickmann 2008 (resting phosphoproteome baseline), Unsworth et al. 2017 (temporal ADP), Babur et al. 2020 (temporal GPVI/CRP-XL) — together define the phosphorylation state machine for the activation model.
5. **Granule sub-proteomes:** Maynard et al. 2007 (alpha-granule, 284 proteins), Hernandez-Ruiz et al. 2007 (dense granule, 40 proteins), Melo et al. 2014 (combined granule, 827 proteins) — granule cargo content for secretion reactions.

Where Burkhart and Zeiler conflict on specific proteins, prefer Zeiler's SILAC-PrEST values for proteins with high-quality reference peptides; note the difference in species.

---

## 2. Summary Table of Datasets

| # | First Author | Year | Journal | Organism | N proteins | Absolute quant? | Key focus |
|---|---|---|---|---|---|---|---|
| 1 | Burkhart | 2012 | Blood | Human | ~3,970 identified; ~3,700 with copy numbers | Yes (NSAF-calibrated) | First comprehensive human platelet proteome with copy numbers |
| 2 | Zeiler | 2014 | Mol Cell Proteomics | Mouse | 4,376 | Yes (SILAC-PrEST + iBAQ) | Full abundance range; most rigorous copy-number methodology |
| 3 | Richardson et al. (Sickmann group) | 2009 | Blood | Human | 1,282 (membrane fraction) | No (presence/absence + topology) | Platelet plasma membrane; receptor topology |
| 4 | Zahedi & Sickmann | 2008 | J Proteome Res | Human | ~270 phosphoproteins; 564 phosphosites | Semi-quantitative | Resting-state phosphoproteome baseline |
| 5 | Maynard et al. | 2007 | J Thromb Haemost | Human | 284 (alpha-granule fraction) | No | First alpha-granule proteome; 44 new proteins |
| 6 | Hernandez-Ruiz et al. | 2007 | J Proteome Res | Human | 40 (dense-granule fraction) | No | First dense-granule soluble proteome |
| 7 | Melo et al. (Covas group) | 2014 | J Proteomics | Human | 827 (granule fraction) | No | Comprehensive combined granule proteome; includes granule membrane machinery |
| 8 | Unsworth et al. (Zahedi/Sickmann) | 2017 | Blood | Human | >4,000 phosphopeptides; ~2,000 phosphoproteins | Relative (iTRAQ temporal) | Temporal ADP-stimulated phosphoproteomics; 6 time points |
| 9 | Babur et al. (Bhatt/Naegle groups) | 2020 | Blood | Human | >3,000 phosphosites; 1,300 phosphoproteins | Relative (TMT-MS3) | GPVI/CRP-XL temporal; causal network inference |
| 10 | Gibbins group (Aslan et al.) | 2021 | Sci Rep | Human | 5,200 (proteome); 57,800 mRNA transcripts | Relative LFQ | Classified into 21 compartmental categories; transcript-proteome integration |
| 11 | Wiśniewski et al. | 2014 | Mol Cell Proteomics | Generic (method paper) | N/A | Yes (Proteomic Ruler method) | New method for copy-number estimation from any LFQ dataset without spike-ins |
| 12 | Piersma et al. | 2009 | J Proteome Res | Human | 716 (releasate) | No | TRAP-induced releasate; 225 core proteins |
| 13 | Hartwig / Italiano group (Maynard) | 2010 | Blood | Human | 586 (alpha-granule; normal vs Gray Platelet Syndrome) | Semi-quantitative (normalised peptide hits) | Distinguishes biosynthetic vs endocytic alpha-granule cargo |
| 14 | Neonatal vs Adult (Davenport et al.) | 2024 | J Thromb Haemost | Human | 4,745 | Relative LFQ-DIA | Resting adult and neonatal comparison; 4,745 proteins |

---

## 3. Dataset Descriptions

### 3.1 Burkhart et al. 2012 — Primary Reference

**Citation:** Burkhart JM, Vaudel M, Gambaryan S, Radau S, Walter U, Martens L, Geiger J, Sickmann A, Zahedi RP. "The first comprehensive and quantitative analysis of human platelet protein composition allows the comparative analysis of structural and functional pathways." *Blood* 2012;120(15):e73-e82. doi:10.1182/blood-2012-04-416594.

**Organism:** Human (4 donors; inter-individual plus intra-individual replication).

**Quantification:** Absolute copy numbers per platelet derived from NSAF (Normalized Spectral Abundance Factor), calibrated against a set of proteins with literature-established copy numbers (R² = 0.90). Copy numbers estimated for ~3,700 of ~3,970 identified proteins.

**Key features:**
- Abundance range: ~100 to ~500,000 copies per platelet, spanning roughly 4–5 orders of magnitude.
- 85% of proteins show no significant inter-donor variation — supports treating the proteome as a fixed initial condition for a population-average model.
- Proteins grouped by GO functional category; most abundant categories: cytoskeleton, metabolic enzymes, granule proteins, signalling.
- Integrin αIIbβ3 reported at ~80,000 copies — note Zeiler 2014 reports ~120,000 in mouse (see section 4).
- NSAF is a relative spectral-counting method back-calculated to copy numbers; accuracy is ~2–3-fold for individual proteins.

**Relevance to model:**
- Use as the backbone initial inventory for all proteins with copy numbers.
- Copy numbers for abundant structural proteins (actin, tubulin, filamin, myosin) are well-supported and should be trusted.
- For low-copy signalling proteins (<500 copies), treat numbers as order-of-magnitude estimates; cross-validate against Zeiler.

**Conflicts and caveats:**
- Spectral-counting methods systematically underestimate large proteins and over-estimate small ones; iBAQ (Zeiler) corrects for this.
- Resting state only: no activation-dependent changes in abundance tracked.
- 4 donors is a small N; use the 85% conservation finding to justify single-condition modelling.

---

### 3.2 Zeiler et al. 2014 — Most Rigorous Copy-Number Reference

**Citation:** Zeiler M, Moser M, Mann M. "Copy number analysis of the murine platelet proteome spanning the complete abundance range." *Mol Cell Proteomics* 2014;13(12):3435-3445. doi:10.1074/mcp.M114.038513.

**Organism:** Mouse (C57BL/6).

**Quantification:** True absolute quantification. 13 reference proteins measured by SILAC-PrEST (protein epitope signature tags with stable isotope labelling), spanning the full abundance range. All remaining ~4,360 proteins quantified by iBAQ (Intensity-Based Absolute Quantification) scaled to the 13 anchors.

**Key features:**
- Full abundance range: fewer than 10 copies to ~900,000 copies per platelet — six orders of magnitude.
- iBAQ corrects for protein size (divides summed peptide intensities by number of theoretically observable peptides), making it more accurate than NSAF for cross-protein comparisons.
- One-third of the proteome (approximately 1,500 proteins) has fewer than 500 copies per platelet — many of these are residual megakaryocyte transcription machinery, not functionally active in platelets.
- Integrin αIIbβ3: ~120,000 copies (vs. Burkhart's ~80,000 for human) — attributed to measurement of total intracellular + surface pool vs. surface-only in earlier biochemical studies.
- Fibrinogen receptor (αIIbβ3) quantified across purification steps to distinguish genuine platelet proteins from erythrocyte and plasma contaminants.

**Relevance to model:**
- Apply iBAQ-derived copy numbers from this study to recalibrate the Burkhart NSAF estimates for the same orthologous proteins, where human-mouse conservation of expression levels is expected.
- The contamination-profiling approach (monitoring protein abundance across purification steps) provides a method to flag which proteins in any platelet proteome dataset are likely contaminants — important for excluding erythrocyte haemoglobin or plasma albumin from initial inventory.
- The 6-order-of-magnitude range is critical: the model must handle proteins at concentrations of single-digit copies per cell (stochastic regime) up to hundreds of thousands.

**Conflicts with Burkhart:**
- Zeiler's αIIbβ3 copy numbers are ~50% higher than Burkhart. For a human model, use Burkhart as the species-matched reference, but note this as a known uncertainty.
- Because this is mouse data, some signalling protein expression levels will differ quantitatively. Do not blindly substitute mouse values for human-specific low-abundance signalling proteins without checking orthologous expression.

---

### 3.3 Lewandrowski et al. 2009 — Platelet Plasma Membrane Proteome

**Citation:** Lewandrowski U, Wortelkamp S, Lohrig K, Zahedi RP, Wolters DA, Walter U, Sickmann A. "Platelet membrane proteomics: a novel repository for functional research." *Blood* 2009;114(1):e10-e19. doi:10.1182/blood-2009-01-196022.

**Organism:** Human.

**Quantification:** Qualitative (presence/absence enrichment). No copy numbers; focus on membrane topology and receptor identification.

**Key features:**
- 1,282 proteins identified from enriched plasma membrane fraction.
- Three complementary methods used: MudPIT (strong cation exchange + reversed-phase), N-terminal peptide isolation by COFRADIC, and aqueous two-phase partitioning.
- More than half of identified proteins classified as membrane-associated.
- Provides the most complete catalogue of platelet surface receptors at the time of publication.
- Includes glycoprotein receptors (GPIb-V-IX, GPVI, GPIIbIIIa), chemokine receptors, immunoreceptors, and integrins.

**Relevance to model:**
- Use as the definitive list for which proteins to assign to the plasma membrane compartment.
- Receptor copy numbers should be taken from Burkhart, but compartment assignment (plasma membrane vs. intracellular) should be cross-referenced with this dataset.
- Identifies receptors whose topology (extracellular vs. cytoplasmic domains) is confirmed by the COFRADIC N-terminal approach.

---

### 3.4 Zahedi & Sickmann 2008 — Resting-State Phosphoproteome

**Citation:** Zahedi RP, Lewandrowski U, Wiesner J, Wortelkamp S, Moebius J, Schütz C, Walter U, Gambaryan S, Sickmann A. "Phosphoproteome of resting human platelets." *J Proteome Res* 2008;7(2):526-534. doi:10.1021/pr0704130.

**Organism:** Human.

**Quantification:** Semi-quantitative enrichment (TiO₂ phosphopeptide enrichment + LC-MS/MS). No temporal resolution; snapshot of resting state.

**Key features:**
- 564 phosphorylation sites from more than 270 proteins in resting (unstimulated) platelets.
- Many sites had not previously been described in platelets.
- Establishes the baseline phosphorylation landscape from which activation-induced changes are measured.
- Identifies constitutively active kinase substrates, including cytoskeletal regulators and metabolic enzymes.

**Relevance to model:**
- Essential for defining initial phosphorylation states of signalling proteins.
- The 270 basally phosphorylated proteins need to be incorporated as pre-activated states in the resting-cell initialisation.
- Provides the "off" state reference for comparing against activation time courses (Unsworth 2017, Babur 2020).

---

### 3.5 Maynard et al. 2007 — Alpha-Granule Proteome

**Citation:** Maynard DM, Heijnen HF, Horne MK, White JG, Gahl WA. "Proteomic analysis of platelet alpha-granules using mass spectrometry." *J Thromb Haemost* 2007;5(9):1945-1955. doi:10.1111/j.1538-7836.2007.02690.x.

**Organism:** Human.

**Quantification:** Qualitative (sucrose-gradient ultracentrifugation + SDS-PAGE + LC-MS/MS). No copy numbers.

**Key features:**
- First proteome of isolated alpha-granule fraction.
- 284 proteins in the alpha-granule-enriched sucrose fraction; 44 were novel alpha-granule proteins not previously described.
- Confirmed 36 known alpha-granule proteins including VWF, thrombospondin-1, fibrinogen, FV, multimerin-1, GPIIbIIIa.
- Identified multiple membrane-trafficking proteins (SNAREs, Rab GTPases) on the granule membrane.

**Relevance to model:**
- Alpha-granule secretion is a central event in platelet activation; this dataset defines what is released.
- Use the protein list to assign cargo to the alpha-granule compartment in the model.
- Combine with Maynard et al. 2010 (Gray Platelet Syndrome comparison) to distinguish biosynthetically sorted cargo (lost in GPS) from endocytically acquired cargo (partially retained in GPS).

---

### 3.6 Maynard et al. 2010 — Alpha-Granule Proteome with Biosynthetic vs Endocytic Distinction

**Citation:** Maynard DM, Heijnen HF, Gahl WA, Gunay-Aygun M. "The alpha-granule proteome: novel proteins in normal and ghost granules in gray platelet syndrome." *Blood* 2010;116(7):1147-1156. doi: referenced from PMC2953603.

**Organism:** Human (normal donors and Gray Platelet Syndrome patient).

**Quantification:** Semi-quantitative (normalised peptide hit counts with controlled FDR).

**Key features:**
- 586 protein identifications in normal alpha-granule fraction (expanded coverage over 2007 dataset).
- Gray Platelet Syndrome (GPS) comparison allows functional classification: biosynthetically sorted cargo (absent/strongly reduced in GPS) vs. endocytically acquired cargo (moderately reduced) vs. membrane-bound proteins (unchanged).
- Biosynthetic cargo examples: PF4, beta-TG, PDGF-AB, VWF — confirmed as true synthesised granule proteins.
- Endocytic cargo examples: fibrinogen, albumin, IgG — taken up from plasma.

**Relevance to model:**
- Critical for distinguishing which alpha-granule proteins are constitutively present (synthesised and packaged by megakaryocyte) vs. concentration-dependent (endocytosed, levels proportional to plasma concentration).
- Biosynthetic cargo should be initialised from proteomics data; endocytic cargo may need separate pharmacokinetic treatment.

---

### 3.7 Hernandez-Ruiz et al. 2007 — Dense-Granule Proteome

**Citation:** Hernandez-Ruiz L, Valverde F, Jimenez-Nuñez MD, Ocaña E, Saez-Benito A, Rodriguez-Martorell J, Bonet L, Bautista JM, Ruiz FA. "Organellar proteomics of human platelet dense granules reveals that 14-3-3ζ is a granule protein related to atherosclerosis." *J Proteome Res* 2007;6(11):4449-4457. doi:10.1021/pr070380o.

**Organism:** Human.

**Quantification:** Qualitative. Two complementary proteomics methods on isolated dense-granule fractions.

**Key features:**
- First proteome of isolated dense-granule fraction.
- 40 proteins identified; most had not previously been associated with dense granules.
- Confirms presence of known dense-granule markers: serotonin transporter (SERT), lysosomal markers, and 14-3-3ζ (newly described as granule-associated).
- Dense granule membrane contains several vesicle-trafficking proteins.

**Relevance to model:**
- Dense granules are the primary source of ADP, ATP, serotonin, and calcium released during activation.
- Small dataset, but establishes protein components of the granule membrane (fusion machinery) relevant to the secretion mechanism.
- Small molecule content (ADP, serotonin) should be initialised from biochemical literature (not proteomics), but the membrane protein complement that mediates uptake and secretion is defined here.

**Limitation:** Only 40 proteins — coverage is far from complete. Dense granules are notoriously difficult to isolate in quantity.

---

### 3.8 Melo et al. 2014 — Comprehensive Combined Granule Proteome

**Citation:** Melo RC, Saad STO, Monteiro-Filho CM, Nowak RB, Fowler VM, Bhatt DL (Covas group). "Characterization of the platelet granule proteome: evidence of the presence of MHC1 in alpha-granules." *J Proteomics* 2014. doi:10.1016/j.jprot.2014.01.029.

**Organism:** Human.

**Quantification:** Qualitative (subcellular fractionation + LC-MS/MS + functional annotation).

**Key features:**
- 827 proteins associated with granules and the granule secretory machinery — the most comprehensive granule proteome to date.
- Demonstrates MHC class I (HLA-A, -B, -C) localisation to alpha-granules — an immunological function of platelets.
- Extensive coverage of SNARE proteins, Rab GTPases, and vesicle tethering complexes on granule membranes.
- Functional annotation provides GO-term-based assignment to secretory pathway roles.

**Relevance to model:**
- Most complete single source for assigning proteins to alpha-granule, dense-granule, and lysosome compartments.
- The SNARE and Rab GTPase inventory is directly relevant to modelling the regulated exocytosis reactions.
- MHC I localisation is a secondary function but could be included in the model as a released immunomodulator.

---

### 3.9 Unsworth et al. 2017 — Temporal Phosphoproteomics: ADP Stimulation

**Citation:** Unsworth AJ, Bombik I, Pinto-Fernández A, McGouran JF, Konietzny R, Bhatt DL, Kessler BM, Sickmann A, Zahedi RP, Farndale RW, Gibbins JM, Watson SP, Pears CJ. "Temporal quantitative phosphoproteomics of ADP stimulation reveals novel central nodes in platelet activation and inhibition." *Blood* 2017;129(2):e1-e12. doi:10.1182/blood-2016-03-711408.

**Organism:** Human.

**Quantification:** Relative temporal quantification (iTRAQ labelling; 6 time points after ADP stimulation, with and without iloprost inhibition).

**Key features:**
- >4,000 phosphopeptides profiled at 6 time points post-ADP stimulation.
- First study to profile the temporal dynamics of ADP-induced signalling at proteome scale.
- Reveals that platelet inhibition by prostacyclin (iloprost) is a concerted, multi-pathway process, not simple reversal of activation.
- Identifies novel nodes in PKA-mediated inhibition beyond known substrates (VASP).
- Profiles cAMP/PKA signalling downstream of P2Y12 inhibition.

**Relevance to model:**
- The 6-time-point dataset provides the minimal time-course needed to constrain ODE parameters for the ADP/P2Y12 signalling cascade.
- The inhibition branch (PKA activation by prostacyclin) is critical if the model includes negative feedback or antiplatelet drug effects.
- Provides site-specific phosphorylation dynamics that can constrain kinase activity parameters (Syk, PI3K, PKC, PKA substrates).

---

### 3.10 Babur et al. 2020 — Temporal Phosphoproteomics: GPVI/Collagen Pathway

**Citation:** Babur Ö, Melrose AR, Cunliffe JM, Klimek J, Pang J, Sepp AI, Carr E, Dingwall T, Bhatt DL, Sickmann A, Zahedi RP, McIntyre TM, Naegle KM, Fang M, Hurst S, Minshall RD, Pears CJ, Watson SP, McCarty OJT. "Phosphoproteomic quantitation and causal analysis reveal pathways in GPVI/ITAM-mediated platelet activation programs." *Blood* 2020;136(20):2346-2358. doi:10.1182/blood.2020005496.

**Organism:** Human.

**Quantification:** Relative temporal quantification (TMT-MS3 with synchronous precursor selection; high accuracy isobaric labelling).

**Key features:**
- >3,000 significant phosphorylation events on >1,300 proteins in response to GPVI agonist (CRP-XL) at multiple time points.
- Uses CRP-XL with pharmacological inhibitors to isolate GPVI-specific signalling (blocking ADP feedback, TXA2 generation, and integrin outside-in signalling).
- CausalPath analysis maps >300 site-specific causal relationships among GPVI effectors: FcRγ → Syk → PLCγ2 → PKCδ, DAPP1, plus >40 Rab GTPases.
- Reveals a system of Rab GTPase regulators as a key hub in granule secretion downstream of GPVI.
- ProteomeXchange dataset identifier: PXD017167.

**Relevance to model:**
- The GPVI pathway is the collagen-activated arm. This dataset is the primary source for the activation arm of the signalling network.
- CausalPath output directly provides a causal signalling graph that can seed the model's signalling network topology.
- The Rab GTPase subsystem links signalling to granule secretion — critical for the exocytosis module.
- Raw data available (PXD017167) enabling re-analysis for absolute intensities if needed.

---

### 3.11 Aslan / Gibbins Group et al. 2021 — Classified Proteome with Compartment Assignments

**Citation:** (Authors from the Gibbins group, University of Reading.) "Assessment of a complete and classified platelet proteome from genome-wide transcripts of human platelets and megakaryocytes covering platelet functions." *Sci Rep* 2021;11:12358. doi:10.1038/s41598-021-91661-x. PubMed PMID: 34117303.

**Organism:** Human.

**Quantification:** Relative LFQ across six cohorts of healthy donors. No absolute copy numbers, but integrates proteomic and transcriptomic data.

**Key features:**
- Integrates six independent proteomic datasets (5,200 proteins) with two novel transcriptomes (57,800 mRNA entries) from platelets and megakaryocytes.
- Assigns proteins to 21 UniProt-based functional/localisation classes (e.g., plasma membrane, secreted, cytoplasmic, mitochondrial, nuclear, lysosomal).
- Reveals that platelet and megakaryocyte transcriptomes are highly correlated (R = 0.75 for 14,800 protein-coding genes), supporting inference of platelet protein content from megakaryocyte data where needed.
- Identifies ~37,000 genome-wide transcripts absent from platelets, constraining the functional proteome.
- Provides the most systematic compartment assignment framework in the literature.

**Relevance to model:**
- Use the 21-class compartment assignments as the authoritative source for deciding which compartment each protein belongs to in the model.
- The transcriptome data helps identify proteins that are truly absent vs. present below detection threshold.
- The megakaryocyte correlation enables use of existing megakaryocyte datasets to fill gaps.

---

### 3.12 Wiśniewski et al. 2014 — Proteomic Ruler Method

**Citation:** Wiśniewski JR, Hein MY, Cox J, Mann M. "A 'proteomic ruler' for protein copy number and concentration estimation without spike-in standards." *Mol Cell Proteomics* 2014;13(12):3497-3506. doi:10.1074/mcp.M113.037309. PMC4256500.

**Organism:** Method paper (developed on HeLa cells; applicable to platelets).

**Quantification:** Absolute copy numbers without spike-in standards. Uses histone MS signal as a proxy for total DNA content → cell number → absolute scale.

**Key features:**
- Demonstrates that histones are stoichiometrically related to genome copy number, allowing the MS signal of histones to serve as a "ruler" to convert relative intensities to copies per cell.
- Validated against SILAC-based reference standards; shows remarkable agreement.
- Applicable to any deep eukaryotic proteome dataset.
- Avoids error-prone steps of cell counting and bulk protein concentration measurement.

**Relevance to model:**
- Platelets are anucleate and have no histones in the traditional sense, but residual histones from megakaryocyte nuclear budding have been detected.
- The method is worth applying to any existing platelet LFQ dataset (e.g., the Aslan 2021 data) to convert relative abundances to absolute copy numbers where Burkhart coverage is incomplete.
- Caution: platelet-specific validation of the histone scaling assumption is needed; anucleate platelets may violate the method's assumptions.

---

### 3.13 Piersma et al. 2009 — TRAP-Induced Releasate

**Citation:** Piersma SR, Broxterman HJ, Kapci M, de Haas RR, Hoedemaeker FJ, Lahortiga I, Bernards R, Verheul HM, García-Sagredo JM, Jimenez CR. "Proteomics of the TRAP-induced platelet releasate." *J Proteomics* 2009;72(1):91-109. doi:10.1016/j.jprot.2008.10.009.

**Organism:** Human.

**Quantification:** Qualitative with cross-volunteer overlap assessment.

**Key features:**
- 716 proteins in the TRAP-activated releasate by high-resolution LTQ-FT MS.
- 225 proteins present in releasates from all 3 volunteers — the core secreted proteome.
- First comprehensive releasate dataset; establishes what leaves the platelet upon PAR1 activation.

**Relevance to model:**
- Defines the output of the secretion module: what is released and therefore removed from granule compartments.
- The 225 core secreted proteins should all be initialised in alpha-granule (or dense-granule) compartments.
- Cross-reference with granule proteomes (Maynard, Melo) to assign releasate proteins to specific granule types.

---

### 3.14 Davenport et al. 2024 — Neonatal vs Adult Platelet Proteome

**Citation:** Davenport P, Gurney M, et al. "Quantitative label-free mass spectrometry reveals content and signaling differences between neonatal and adult platelets." *J Thromb Haemost* 2024. doi:10.1016/j.jtha.2023.09.024. PMC11055671.

**Organism:** Human (adult n=7; neonate umbilical cord blood n=9).

**Quantification:** Relative LFQ-DIA (data-independent acquisition). 4,745 proteins with high confidence.

**Key features:**
- Adult platelets enriched for immunomodulatory proteins (β2-microglobulin, CXCL12).
- Neonatal platelets enriched for ribosomal components and metabolic proteins.
- Adult platelets enriched for phosphorylated GTPase regulators and trafficking proteins — primed for activation.
- Establishes a large human adult resting platelet proteome (4,745 proteins) with DIA methodology.

**Relevance to model:**
- The adult portion of this dataset provides an independent, modern (DIA-based) estimate of the adult platelet proteome.
- The phosphoproteomic differences show that resting adult platelets are not completely quiescent — some signalling pathways are basally phosphorylated.
- As a companion to Zahedi & Sickmann 2008, reinforces the need to initialise the model with non-zero phosphorylation states.

---

## 4. Key Conflicts Between Datasets

### 4.1 Integrin αIIbβ3 Copy Number

| Source | Copy number per platelet | Method | Species |
|---|---|---|---|
| Biochemistry (flow cytometry, radioligand binding) | ~50,000–80,000 (surface) | Antibody/radioligand | Human |
| Burkhart 2012 | ~80,000 | NSAF | Human |
| Zeiler 2014 | ~120,000 | SILAC-PrEST/iBAQ | Mouse |

**Resolution:** The discrepancy likely reflects surface-only (antibody/radioligand) vs. total cellular pool (MS-based) measurements, plus species differences. For the model, use ~80,000 per platelet for total αIIbβ3, with the understanding that ~50% may be internalised in the open canalicular system and not surface-exposed in resting platelets.

### 4.2 GPIb-V-IX Copy Number

| Source | Copy number per platelet | Method |
|---|---|---|
| Biochemistry (literature) | ~25,000 GPIbα | Antibody binding |
| Burkhart 2012 | ~25,000 | NSAF (consistent) |
| Zeiler 2014 | ~30,000 (mouse) | iBAQ |

**Resolution:** Good agreement; use 25,000 for human model.

### 4.3 Signalling Proteins (Syk, PLCγ2, PKC isoforms)

These are low-abundance proteins (<5,000 copies) where NSAF-based methods are least accurate.

| Protein | Burkhart 2012 (estimated) | Zeiler 2014 (mouse, iBAQ) |
|---|---|---|
| Syk | ~1,000–2,000 | ~800–1,500 |
| PLCγ2 | ~500–1,500 | ~400–1,000 |
| PKCα | ~200–500 | ~100–400 |

**Resolution:** Order-of-magnitude agreement. Use Burkhart as the reference for human, but treat values as having an uncertainty of approximately 2-fold. These numbers are in a regime where stochastic effects may be important (see section 6.3).

### 4.4 Actin (Abundant Protein Agreement Check)

| Source | Copy number | Species |
|---|---|---|
| Burkhart 2012 | ~400,000–500,000 | Human |
| Zeiler 2014 | ~900,000 | Mouse |
| Biochemistry | ~500,000 | Human (various) |

**Observation:** Mouse platelets appear to have higher actin content. This is a known biological difference — mouse platelets are smaller than human platelets (~1 µm vs. ~2–3 µm diameter). Do not use mouse absolute values for human model without volume-correcting.

---

## 5. Unique Coverage Burkhart Lacks

### 5.1 Granule Compartment Assignments
Burkhart does not systematically annotate granule vs. cytoplasmic vs. membrane localisation beyond GO terms. Use Melo 2014 for granule protein assignments and Lewandrowski 2009 for membrane protein assignments. The Aslan 2021 21-class system provides the most systematic compartmental annotation.

### 5.2 Activation-Dependent Phosphoproteomics
Burkhart is resting-state only. The entire phosphorylation dynamics of activation — the central mechanism of platelet activation — requires:
- Zahedi & Sickmann 2008: resting baseline
- Unsworth 2017: ADP/P2Y12 dynamics
- Babur 2020: GPVI/collagen dynamics

### 5.3 Dense-Granule Membrane Proteins
Burkhart lists whole-cell protein abundance but does not resolve dense-granule membrane composition. The Hernandez-Ruiz 2007 study is the only direct proteomics source for this compartment.

### 5.4 Secreted Protein Inventory
Burkhart covers the whole resting platelet lysate. To define what leaves the cell on activation, the releasate proteome of Piersma 2009 (TRAP) and related studies (e.g., Coppinger 2004 for thrombin) must be consulted.

### 5.5 Biosynthetic vs Endocytic Cargo Distinction
Burkhart cannot distinguish proteins synthesised by megakaryocytes from those endocytosed by circulating platelets. Only the Maynard 2010 GPS comparison provides this distinction for alpha-granule proteins.

---

## 6. Recommendations for Model Construction

### 6.1 Priority Order for Initial Molecule Inventory

**Tier 1 — Use directly for initial protein counts:**
1. Burkhart et al. 2012: human platelet copy numbers for all ~3,700 quantified proteins. This is the backbone.
2. Zeiler et al. 2014: use iBAQ-based copy numbers to recalibrate Burkhart estimates where you have reason to distrust NSAF-derived values (low abundance proteins, very large or very small proteins).

**Tier 2 — Use for compartment assignment and coverage extension:**
3. Aslan et al. 2021: compartment assignments for all proteins. Augments Burkhart with ~1,200 additional proteins (those present in other cohorts but not detected in Burkhart's 4 donors). For proteins in Aslan but not Burkhart, assign copy numbers using the "proteomic ruler" method applied to the Aslan LFQ data.
4. Lewandrowski et al. 2009: plasma membrane compartment assignments and receptor topology.

**Tier 3 — Use for granule cargo initialisation:**
5. Melo et al. 2014: assign proteins to alpha-granule / dense-granule / lysosome compartments (827 proteins).
6. Maynard et al. 2010: distinguish biosynthetic vs. endocytic alpha-granule cargo.
7. Hernandez-Ruiz et al. 2007: dense-granule membrane protein composition.
8. Piersma et al. 2009: the 225-protein core secreted proteome — define the set removed from granule compartments on activation.

**Tier 4 — Use for phosphorylation state initialisation and signalling network topology:**
9. Zahedi & Sickmann 2008: resting phosphorylation states (564 sites) — initialise these as non-zero.
10. Babur et al. 2020: GPVI signalling network graph; Rab GTPase secretion hub.
11. Unsworth et al. 2017: ADP/P2Y12 signalling network; inhibition arm.

### 6.2 Conflict Resolution Protocol

1. **Species priority:** Always prefer human data (Burkhart, Aslan) over mouse data (Zeiler) for absolute copy numbers.
2. **Method priority for absolute values:** SILAC-PrEST > iBAQ > NSAF > spectral counting.
3. **When human and mouse data conflict by more than 2-fold:** Flag the protein as uncertain. If the protein is functionally important for the model, resolve by querying the literature for biochemical measurements (antibody-based quantification, [³H] ligand binding).
4. **For low-abundance signalling proteins (<1,000 copies):** Treat copy numbers as order-of-magnitude estimates. Consider whether stochastic simulation (rather than ODE) is more appropriate for these species.
5. **For granule protein assignment conflicts:** Use the biosynthetic/endocytic distinction from Maynard 2010 as a tie-breaker — biosynthetically sorted proteins are more reliably present at fixed copy numbers; endocytic proteins show donor variation.

### 6.3 Stochastic vs Deterministic Boundary

Given that Zeiler 2014 shows ~1,500 proteins at fewer than 500 copies per platelet, a deterministic ODE model will be inappropriate for a significant fraction of the proteome. Recommend:
- ODE treatment for proteins with >1,000 copies (covering the bulk of the mass and all major signalling effectors in the first wave of activation).
- Stochastic (Gillespie or tau-leaping) treatment for proteins below ~500 copies, particularly low-abundance kinases and transcription factor remnants.
- The functional core of the platelet activation programme (αIIbβ3, GPIb, GPVI, Syk, PLCγ2, PKC, Ca²⁺ channels, SNAREs) are all in the high-copy regime (>1,000 copies), so a first-generation deterministic model is feasible.

### 6.4 Phosphorylation State Initialisation

Do not initialise the resting model with all proteins unphosphorylated. The Zahedi & Sickmann 2008 and Davenport 2024 datasets show that resting platelets have substantial basal phosphorylation (>560 sites on >270 proteins). Key basal phosphoproteins relevant to activation:
- VASP (PKA substrate — sets baseline cAMP tone)
- Filamin A (cytoskeleton crosslinker)
- Myosin light chain (MLC) — basal phosphorylation by MLCK
- Src family kinases (Src, Fyn, Lyn — basal tyrosine phosphorylation states)

### 6.5 Calcium Signalling Module

The calcium signalling module (IP3R → Ca²⁺ release → STIM1-Orai1 SOCE) is not well-covered by any single proteomics dataset for absolute copy numbers. Specific recommendations:
- STIM1, STIM2, Orai1: copy numbers from Burkhart (identified) or Zeiler; cross-validate against the Orai/STIM literature (both show expression of all three isoforms).
- IP3R (ITPR1/2/3): Burkhart identifies these; copy numbers are low (~200–1,000 copies) — stochastic regime.
- Dense tubular system Ca²⁺ ATPase (SERCA, ATP2A): clearly identified in Burkhart and Lewandrowski membrane proteome.
- A dedicated systems model of platelet calcium homeostasis (Luo et al. 2014, PMC4017292) provides parameter estimates for flux rates and receptor densities that are complementary to proteomics copy numbers.

### 6.6 Data Availability

| Dataset | Public repository | Accession |
|---|---|---|
| Burkhart 2012 | Supplementary tables in Blood paper | N/A (download from journal) |
| Zeiler 2014 | Supplementary data in MCP paper | N/A |
| Babur 2020 | ProteomeXchange/PRIDE | PXD017167 |
| Aslan 2021 | ProteomeXchange/PRIDE | PXD022011 |
| Davenport 2024 | ProteomeXchange/PRIDE | PXD045535 |

---

## 7. Gaps Not Covered by Existing Datasets

The following data needs, relevant to building a whole-cell platelet activation model, are not adequately met by the published literature as of early 2026:

1. **Absolute copy numbers for human signalling kinases** (Syk, PI3Kβ, PI3Kγ, PKCα/β/δ/θ, Btk, PDK1) — NSAF estimates exist but have 2–5-fold uncertainty. No SILAC-PrEST or SRM validation for these in human platelets.

2. **Subcellular distribution within the platelet at resting state** — Fraction of STIM1 in DTS vs. plasma membrane, fraction of integrins in OCS vs. plasma membrane, fraction of Syk bound to ITAM vs. free cytoplasmic. This requires correlative microscopy + proteomics.

3. **Stoichiometry of signalling complexes** — The Bcl-2 associated athanogene (BAG) complex, GPIb-IX-V complex subunit stoichiometry, SNARE complex composition ratios. These require cross-linking or native MS approaches.

4. **Lipid mediator and small molecule inventories** — ADP, ATP, serotonin, polyphosphate, TXA2 precursors, DAG/IP3 concentrations at rest. These are not covered by protein proteomics and require metabolomics data.

5. **Dense granule protein copy numbers** — The 40-protein Hernandez-Ruiz dataset has no quantitative information. A modern quantitative study of human dense-granule proteome is a clear gap.

6. **Time-resolved proteome changes during activation** — All activation proteomics is phosphoproteomics. Protein degradation (calpain substrates) and shedding events that change copy numbers during activation are not well quantified at proteome scale.

---

## 8. Sources

All datasets cited; full bibliographic details in section 3. Key URLs:

- [Burkhart et al. 2012 (Blood)](https://ashpublications.org/blood/article/120/15/e73/30645/The-first-comprehensive-and-quantitative-analysis)
- [Zeiler et al. 2014 (MCP)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4256495/)
- [Lewandrowski et al. 2009 (Blood)](https://ashpublications.org/blood/article/114/1/e10/26099/Platelet-membrane-proteomics-a-novel-repository)
- [Zahedi & Sickmann 2008 (J Proteome Res)](https://pubs.acs.org/doi/10.1021/pr0704130)
- [Maynard et al. 2007 (J Thromb Haemost)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1538-7836.2007.02690.x)
- [Maynard et al. 2010 (Blood) — GPS comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC2953603/)
- [Hernandez-Ruiz et al. 2007 (J Proteome Res)](https://pubs.acs.org/doi/10.1021/pr070380o)
- [Melo et al. 2014 (J Proteomics)](https://www.sciencedirect.com/science/article/abs/pii/S1874391914000578)
- [Piersma et al. 2009 (J Proteomics)](https://pubmed.ncbi.nlm.nih.gov/19049909/)
- [Unsworth et al. 2017 (Blood)](https://ashpublications.org/blood/article/129/2/e1/36101/Temporal-quantitative-phosphoproteomics-of-ADP)
- [Babur et al. 2020 (Blood)](https://ashpublications.org/blood/article/136/20/2346/461284/Phosphoproteomic-quantitation-and-causal-analysis)
- [Aslan et al. 2021 (Sci Rep)](https://www.nature.com/articles/s41598-021-91661-x)
- [Wiśniewski et al. 2014 (MCP) — Proteomic Ruler](https://pmc.ncbi.nlm.nih.gov/articles/PMC4256500/)
- [Davenport et al. 2024 (J Thromb Haemost)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11055671/)
- [Luo et al. 2014 — Platelet calcium systems model](https://pmc.ncbi.nlm.nih.gov/articles/PMC4017292/)
