# Honoré and Vorum (2000) — The CREC family of low-affinity ER/secretory Ca²⁺-binding proteins

- **Full citation:** Honoré B, Vorum H. "The CREC family, a novel family of multiple EF-hand, low-affinity Ca²⁺-binding proteins localised to the secretory pathway of mammalian cells." *FEBS Letters* 466:11–18, 2000.
- **Type:** review (minireview)
- **Local PDF:** source-info/calcium-papers/Honoré and Vorum - 2000 - The CREC family, a novel family of multiple EF-hand, low-affinity Ca2+-binding proteins localised to.pdf
- **Used in the model for:** Ca²⁺-binding parameters (EF-hand count, low affinity, localisation) of the CREC-family buffers calumenin (CALU), reticulocalbin (RCN1) and reticulocalbin-2 (RCN2) in the multi-buffer DTS store (issue #25, Phase 3).

## What the paper does
Reviews the CREC family — mammalian members **reticulocalbin (RCN1), ERC-55/TCBP-49/E6BP (RCN2), Cab45 (SDF4), calumenin (CALU) and crocalbin/CBP-50** — multiple-EF-hand, low-affinity Ca²⁺-binding proteins of the secretory pathway. Covers structure/localisation, the EF-hand motif and its Ca²⁺ affinity, gene structure/evolution, and disease associations (cancer, snake-venom toxin mediation, amyloid).

## Key findings / values the model uses
- Mammalian CREC proteins are **acidic pro-proteins (315–362 aa, ~37–42 kDa; mature ~35–40 kDa)** with N-terminal signal sequences and **six, or seven (calumenin, crocalbin), EF-hand domains**.
- **Ca²⁺ affinity is characteristically low: dissociation constants ~10⁻⁴–10⁻³ M (up to the mM range)** — unusual for EF-hand proteins.
- For **human calumenin**, the one solution study (ref 13 = Vorum 1998) found **no high-affinity site**; all seven EF-hands bind Ca²⁺ with similar low affinity, **Kd ≈ 0.6×10⁻³ M (0.6 mM) at 37 °C**. RCN1 and ERC-55 likewise lack a high-affinity site (authors' unpublished).
- Localisation: reticulocalbin and ERC-55 strictly ER (C-terminal HDEL retrieval); Cab45 strictly Golgi (HEEF); **calumenin** distributed through the secretory pathway (HDEF, a weaker retrieval signal) and **uniquely secreted**.
- CREC proteins are **low-abundance** (ERC-55 below Coomassie/silver detection in HeLa, far below endoplasmin/BiP/PDI/calreticulin), so a storage/buffer role is "not very likely" to dominate; their low-affinity binding fits Ca²⁺-regulated roles in the high-Ca²⁺ secretory lumen/extracellular space.

## Relevance to the platelet WCM
Supplies the buffer parameters for the CALU/RCN1/RCN2 species in the planned multi-buffer DTS store: ~6–7 low-affinity EF-hand sites per molecule, Kd in the 0.1–1 mM range. The "low abundance" caveat suggests these contribute less to the DTS Ca²⁺ budget than the high-capacity BiP/CRT/GRP94 chaperones — relevant when apportioning the store.

## Caveats / notes
A review; the only directly-measured affinity it reports is calumenin's (from Vorum 1998). Mammalian CREC localisation is ER/secretory-pathway generic, not platelet-DTS-specific. Use platelet proteome copy numbers (Burkhart/Huang) for the in-platelet quantities, which are not given here. Cross-referenced in `reports/data/calcium-data-provenance.md` (CREC review; CALU/RCN1/RCN2; Kd up to mM).
