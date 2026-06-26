# Burkhart et al. (2012) — Quantitative human platelet proteome with per-platelet copy numbers

- **Full citation:** Burkhart JM, Vaudel M, Gambaryan S, Radau S, Walter U, Martens L, Geiger J, Sickmann A, Zahedi RP. "The first comprehensive and quantitative analysis of human platelet protein composition allows the comparative analysis of structural and functional pathways." *Blood* 120(15):e73–e82, 2012.
- **Type:** proteome (quantitative mass spectrometry)
- **Local PDF:** source-info/calcium-papers/Burkhart et al. - 2012 - The first comprehensive and quantitative analysis of human platelet protein composition allows the c.pdf
- **Used in the model for:** primary source of per-platelet copy numbers for the species inventory (`reports/params/species-v0.6.tsv`), especially the calcium-pathway channels, pumps and receptors.

## What the paper does
Highly purified platelets from 4 healthy donors (<1 leukocyte per 10⁶ platelets) were analysed by quantitative MS. The authors identified ~4000 unique proteins and estimated per-platelet copy numbers for ~3700 of them, using NSAF (normalised spectral abundance factor) calibrated against 24 reference proteins with literature copy numbers (Table 1). They estimate ~20 million protein molecules per platelet (~1.5 mg protein / 10⁹ platelets) and ~80–85% proteome coverage.

## Key findings / values the model uses
Copy numbers stated directly in the calcium-signalling paragraph (p. e78) and Table 1:
- **IP₃ receptors:** ITPR1 2400, ITPR2 1700, ITPR3 750 copies (IRAG 3500).
- **PLC isoforms:** PLCB2 2500, PLCB3 1700, PLCB4 1000, PLCG2 2000 copies.
- **SOCE:** STIM1 7400, CRACM1/Orai1 1700, TRPC6 1100 copies.
- **Ca²⁺ pumps/transporters:** PMCA4 640, ATPase 2C1 (SPCA1) 2200, Na⁺/Ca²⁺ exchanger SLC8A3 580; mitochondrial MCU 5900, MICU1 1400; ER/SR Ca²⁺-ATPases SERCA2 9000 and SERCA3 16,300.
- **Ionotropic:** P2X1 1400 (the only ligand-gated Ca²⁺ channel found).
- **Integrins (Table 1):** ITGA2B (integrin αIIb / CD41) 83,300, ITB3 (integrin β3 / CD61) 64,200 — "almost equally expressed."
- Reference set (Table 1): actin ~2.19×10⁶, fibrinogen β 88,900, P-selectin 8900, PKCβ 9700, etc.
- Copy numbers correlate with literature reference values at R² ≈ 0.90; spectral-counting estimates can err up to ~200% in individual cases (membrane proteins, heavily-modified proteins).

## Relevance to the platelet WCM
This is THE provenance source for the channel/pump/receptor copy numbers seeded into `species-v0.6.tsv` and used as initial bulk-molecule counts. The αIIb/β3 figures underpin the integrin-activation module's ~80,000 assembled-heterodimer pool. The very high SERCA:IP3R ratio it reports supports the store-limited behaviour the calcium core exhibits.

## Caveats / notes
Copy numbers are estimates, not absolute; membrane and multi-TMD proteins (e.g. P2X1, GPCRs) are systematically under-counted, so model values should not be treated as exact. Many GPCRs were only detected in the authors' earlier membrane-proteome study (P2Y1 ranked ~482nd; P2Y12 ~296th — copy number not given in the text). Cross-referenced in `reports/data/calcium-data-provenance.md` (Zotero itemKey 3KX7BHEV, verified 2026-04-23).
