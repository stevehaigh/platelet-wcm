# Huang et al. (2021) — Classified, transcriptome-extended human platelet proteome

- **Full citation:** Huang J, Swieringa F, Solari FA, Provenzale I, Grassi L, De Simone I, Baaten CCFMJ, Cavill R, Sickmann A, Frontini M, Heemskerk JWM. "Assessment of a complete and classified platelet proteome from genome-wide transcripts of human platelets and megakaryocytes covering platelet functions." *Scientific Reports* 11:12358, 2021. https://doi.org/10.1038/s41598-021-91661-x
- **Type:** proteome (meta-analysis / integration of proteomes + genome-wide transcriptomes)
- **Local PDF:** source-info/calcium-papers/Huang et al. - 2021 - Assessment of a complete and classified platelet proteome from genome-wide transcripts of human plat.pdf
- **Used in the model for:** secondary / corroborating source for the species inventory; classifies platelet proteins by function and confirms presence of pathway components.

## What the paper does
The authors integrated established platelet proteomes from six cohorts (22 healthy subjects, same analytical workflow) into a merged proteome of 5211 identified proteins, then combined this with genome-wide platelet and megakaryocyte transcriptomes (~57.8k mRNAs) from the Blueprint consortium. Proteins with relevant transcript levels (log2fpkm ≥ 0.20) were assigned to 21 UniProt-based function classes by intracellular localisation/function. For 3629 proteins, copy numbers (carried over from the cohort proteomes) were also available.

## Key findings / values the model uses
- **Identified proteome** set at **5050 proteins** (with relevant mRNA); a validation cohort of pooled platelets from 30 subjects added 954 previously unidentified proteins (Fig. 2).
- A **prediction model** projects a total achievable platelet proteome of **~10,000 proteins**.
- Protein vs transcript level correlation is weak (R = 0.25, Fig. 5A,B): a "triangular" pattern showing transcript level restricts the *maximum* attainable copy number but does not determine it. Low-abundance proteins (<500 copies/platelet) cluster at the floor.
- Three restraining factors for non-detection: low copy number (43% identified), low mRNA >20% (45%), retention in the megakaryocyte (20%); baseline ~65% for other classes.
- **Hemostasis/thrombosis relevance:** a Reactome-based network incorporated 1.3k identified proteins (median ~2200 copies); 138 proteins linked to platelet-related disorders (Fig. 8); 124 of those disorder genes had relatively high platelet copy numbers (median ~22.8k).
- This paper carries copy numbers forward from the cohort proteomes (Burkhart 2012 is cohort 1) rather than re-measuring them; specific Ca²⁺-channel numbers are not re-quoted in the main text.

## Relevance to the platelet WCM
Confirms that the channels, pumps, receptors and integrin subunits used in the model are genuinely part of the classified platelet proteome and assigns them to function classes (e.g. C10 membrane receptors & channels). Useful as an independent corroboration layer above Burkhart for "is this protein really present in platelets" decisions and for the disorder/Glanzmann framing.

## Caveats / notes
This is an integrative/classification study; the per-platelet copy numbers it reports are inherited from the source cohorts (notably Burkhart 2012), so Burkhart remains the authoritative numeric source for the inventory. The weak protein–transcript correlation it documents is a caution against ever inferring platelet copy numbers from mRNA data. Pages 9–11 (later discussion/references) were not read; summary is grounded in pages 1–8 read from source.
