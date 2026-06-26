# Argon and Simen (1999) — GRP94 / HSP90B1, an abundant ER chaperone and Ca²⁺-buffer

- **Full citation:** Argon Y, Simen BB. "GRP94, an ER chaperone with protein and peptide binding properties." *Seminars in Cell & Developmental Biology* 10:495–505, 1999.
- **Type:** review
- **Local PDF:** source-info/calcium-papers/Argon and Simen - 1999 - GRP94, an ER chaperone with protein and peptide binding properties.pdf
- **Used in the model for:** Ca²⁺-binding stoichiometry / affinity for the GRP94 (HSP90B1) component of the multi-buffer DTS luminal Ca²⁺ store (issue #25, Phase 3).

## What the paper does
Reviews GRP94 (also endoplasmin / CaBP4 / gp96 / ERp99), the ER-resident member of the HSP90 family. Covers its gene/expression (stress-induced in parallel with BiP), structure (mature murine protein 782 residues, homodimer, C-terminal KDEL retention), post-translational modification, drug/ATP binding, restricted set of protein substrates (late folding intermediates), and its peptide-binding / antigen-presentation (tumour vaccine) role.

## Key findings / values the model uses
- GRP94 is **one of the most abundant ER proteins, ~5–10% of luminal content, estimated concentration ~10 mg/mL** (Introduction).
- It is a **low-affinity, high-capacity calcium-binding protein** and "one of the important ER calcium buffer proteins" (Calcium binding section).
- **15 calcium-binding sites per molecule: 4 of moderate affinity (Kd ≈ 2 µM) + 11 of low affinity (Kd ≈ 600 µM).**
- Calcium binding is attributed to highly negatively charged regions distributed through the sequence (no obvious EF-hand motifs).
- 100 nM free Ca²⁺ induces a conformational change (helical content drops from ~40% to ~34%); Ca²⁺ may regulate GRP94's in vivo protein interactions.

## Relevance to the platelet WCM
Supplies the per-molecule Ca²⁺-binding parameters (15 sites; 4 × Kd ~2 µM, 11 × Kd ~600 µM) for representing GRP94/HSP90B1 as one species in the planned multi-buffer DTS store. Combined with its high abundance, this makes GRP94 a meaningful luminal Ca²⁺ sink alongside BiP and the CREC proteins.

## Caveats / notes
This is a 1999 review; the stoichiometry/affinity figures are cited to its ref 3 (Van et al. 1989, J Biol Chem 264:17494). The platelet-specific abundance of GRP94 (vs the ~10 mg/mL generic ER figure) is not given here — use Burkhart/Huang proteome copy numbers for the platelet quantity. Cross-referenced in `reports/data/calcium-data-provenance.md` ("15 sites... 4 medium-affinity Kd ~2 µM + 11 low-affinity Kd ~600 µM").
