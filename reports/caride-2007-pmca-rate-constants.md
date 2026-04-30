# Caride 2007 PMCA + CaM rate constants for Phase 1

**Source:** Caride et al. 2007, J Biol Chem 282:25640–25648 (Zotero: YLIICDWW). Pulled via Zotero MCP since the PDF isn't in `source-info/calcium-papers/`.

**Isoform of interest:** PMCA4b (the better-studied isoform; same isoform Dolan 2014 cites). Burkhart 2012 reports total PMCA4 = 640 copies/platelet, no isoform breakdown.

---

## 1. The Caride simplified Ca²⁺-signal model (Table 3)

Caride simulates a hypothetical cell's Ca²⁺ spike by combining:

- ER release on agonist binding (steps 1–3, fitted)
- Plasma-membrane Ca²⁺ entry (step 3, fitted as "store-operated")
- CaM-independent PMCA activity (steps 4–5)
- Ca²⁺ binding to CaM (steps 6–7)
- CaM-mediated PMCA activation (steps 8–11)
- SERCA refilling the store (steps 12–14)

The "kinetic model for CaM activation of PMCA4b" is the part we want — steps 4–11.

---

## 2. Rate constants for PMCA4b (the platelet-relevant set)

All values from Caride 2007 Table 3, "Fit PMCA4b" column:

| Step | Reaction | Rate constant | Value |
|------|----------|---------------|-------|
| 4 | PMCA + Ca²⁺_cyt ⇌ PMCA·Ca | k4 (µM⁻¹·s⁻¹) | 10 |
| 4 | reverse | k4r (s⁻¹) | 50 |
| 5 | PMCA·Ca → PMCA + Ca²⁺_ex | k5 (s⁻¹) | **5.5**  ← basal turnover (currently used) |
| 6 | CaM + 2 Ca²⁺ ⇌ Ca₂·CaM | k6 (µM⁻²·s⁻¹) | 2.669 |
| 6 | reverse | k6r (s⁻¹) | 2.682 |
| 7 | Ca₂·CaM + 2 Ca²⁺ ⇌ Ca₄·CaM | k7 (µM⁻²·s⁻¹) | 170.4 |
| 7 | reverse | k7r (s⁻¹) | 1.551 |
| 8 | PMCA + Ca₄·CaM ⇌ Ca₄·CaM·PMCA | k8 (µM⁻¹·s⁻¹) | 0.2 |
| 8 | reverse | k8r (s⁻¹) | 0.0008 |
| 9 | Ca₄·CaM·PMCA + Ca²⁺_cyt ⇌ Ca₄·CaM·PMCA·Ca | k9 (µM⁻¹·s⁻¹) | 50 |
| 9 | reverse | k9r (s⁻¹) | 10 |
| 10 | Ca₄·CaM·PMCA·Ca → Ca₄·CaM·PMCA + Ca²⁺_ex | k10 (s⁻¹) | **30**  ← CaM-activated turnover (~5.5× basal) |
| 11 | Ca₄·CaM·PMCA ⇌ PMCA·CaM + 4 Ca²⁺ | k11 (s⁻¹) | 10 |
| 11 | reverse | k11r (µM⁻⁴·s⁻¹) | 0.0007332 |

For comparison, PMCA4a values where they differ:
- k5 = 12 s⁻¹ (vs 5.5 for 4b)
- k8 = 0.8 µM⁻¹·s⁻¹ (4× faster CaM association)
- k8r = 0.02 s⁻¹ (25× faster CaM dissociation)
- k11 = 6.2 s⁻¹

PMCA4b is what we want.

---

## 3. State diagram for implementation

The 5-state PMCA cycle plus 3 CaM states:

```
         PMCA states:              CaM states:
                                     CaM_free  ←──────┐
              (4)                        ↓ (6)        │
   PMCA  ←────────→  PMCA·Ca           Ca₂·CaM        │
     ↑(8)              (5) → Ca²⁺_ex     ↓ (7)        │ (8 reverse)
     ↓                                  Ca₄·CaM ──────┘
   Ca₄·CaM·PMCA                         (active)
     ↓ (9)                                ↑
   Ca₄·CaM·PMCA·Ca                       │
     (10) → Ca²⁺_ex                      │
     ↑(11)                                │
   PMCA·CaM ──── (steps not detailed)────┘
```

This adds 6 new sub-states beyond the current 2 (PMCA, PMCA·Ca):
- `Ca4_CaM_PMCA[pl]`        (CaM-activated, empty)
- `Ca4_CaM_PMCA_Ca[pl]`     (CaM-activated, Ca²⁺-loaded)
- `PMCA_CaM[pl]`            (CaM bound, no Ca²⁺ on CaM)
- Plus 3 CaM-only states: `CaM_free[c]`, `Ca2_CaM[c]`, `Ca4_CaM[c]`

(`PMCA·CaM` is functionally similar to a "deactivated" CaM-bound form; can keep or skip depending on whether step 11 matters at our timescales.)

---

## 4. Initial counts from Dolan 2014 Table S1

| Sub-state | Dolan IC count | Notes |
|-----------|----------------|-------|
| PMCA (free) | 765 | Already in model |
| Ca²⁺·PMCA | 4 | Already in model |
| PMCA(Ca₄·CaM) | ~0 | New |
| PMCA(Ca₄·CaM)·Ca²⁺ | ~0 | New |
| PMCA·CaM | ~0 | New |
| CaM (free) | 20,465 | New |
| Ca²⁺₂·CaM | 15 | New |
| Ca²⁺₄·CaM | 1 | New |
| (CaM bound to PMCA) | ~0 | New |

Total CaM = 20,481 (Dolan Table S1, verified).
Total PMCA = 769 (Dolan), of which 765+4 = 769 are unbound to CaM at rest.

---

## 5. Caveats / things I had to interpret

1. **Step 11 is not fully described.** Caride lists `k11=10 s⁻¹` and `k11r=0.0007 µM⁻⁴·s⁻¹` — the µM⁻⁴ implies a 4th-order Ca²⁺ rebinding. Likely scheme: `Ca₄·CaM·PMCA ⇌ PMCA·CaM + 4 Ca²⁺` (slow Ca²⁺ release from CaM while CaM stays bound to PMCA). The original kinetic detail is in Penheiter 2003 (cited; not in Zotero). For Phase 1 we can leave step 11 out — it's slow (10 s⁻¹) and only matters during long-timescale relaxation.

2. **CaM Ca²⁺ binding (steps 6–7) is two-lobe, not Hill.** Caride splits CaM Ca²⁺ binding into two cooperative pairs (N-lobe and C-lobe), each binding 2 Ca²⁺. This is more accurate than a single 4-Ca cooperative step. The k6 (2.669 µM⁻²·s⁻¹) is the slow lobe; k7 (170.4 µM⁻²·s⁻¹) is the fast lobe.

3. **k4 / k4r / k5 are the same numbers we already use** for basal PMCA. So the "basal-only PMCA" path doesn't change; we add the CaM-activated path as an additional cycle.

4. **k8r = 0.0008 s⁻¹ is very slow** (CaM dissociation from PMCA). Once Ca₄·CaM binds PMCA, it stays bound for ~20 minutes. This is what makes PMCA "memory" of recent Ca²⁺ spikes physiologically meaningful.

5. **CaM concentration in the platelet** = 20,481 / (6 fL × N_A) ≈ 5.66 µM. Caride's CHO cells have ~1.94 µM CaM; platelets have 3× more.

---

## 6. Implementation notes for Phase 1b–1d

- All Caride rates assume concentrations in **µM**. Conversion to count form uses our existing `_UM_PER_COUNT_CYT = 2.77×10⁻⁴ µM/count`.
- `k8` (PMCA + CaM bimolecular) is a count⁻¹·s⁻¹ rate when applied to two count quantities. `k8 [µM⁻¹·s⁻¹] × _UM_PER_COUNT_CYT [µM/count] = 2.77×10⁻⁵ count⁻¹·s⁻¹` for the PMCA/CaM bimolecular step. Both species live in cytosolic volume, so this is consistent.
- For the CaM Ca²⁺-binding (steps 6, 7), the rate `k6 × CaM_count × ca_cyt²` has units of `µM⁻²·s⁻¹ × count × µM²` = count/s ✓.

Phase 1d will need `K_CAM_PMCA = {...}` and `K_CAM = {...}` dicts mirroring the existing `K_SOCE` style.

---

*Phase 1a complete (2026-04-30). Phase 1b → 1d follow once user confirms.*
