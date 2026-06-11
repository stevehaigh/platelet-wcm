---
title: "Lab book — 2026-06-10: out/ tidy-up, decay-stub regeneration, PMCA/MCU perturbation figures (#53)"
---

# Lab book — 2026-06-10: data tidy-up + the two mechanistic-finding figures

## Context — bridging the gap since 2026-05-22

The lab-book trail went cold after the 2026-05-22 dose-sweep entry. The
intervening work was infrastructure and writing, not new biology:

- **#32 kinetics-as-data** landed — ~150 rate constants carved to
  `reports/params/calcium-v0.5.toml`, the 63-species inventory to
  `species-v0.5.tsv`, both loaded at import. Model relabelled **v0.5**, but
  the calcium biology is unchanged from **v0.4.1** (2026-05-12).
- The **MRes thesis** scaffold was built in `reports/thesis/`: `draft-thesis.qmd`,
  `reflective_document.qmd`, four supplement skeletons, the COPASI/VCell
  supplement, and the advisory `review-part-2.md` (assessment split 70/20/10,
  marking-descriptor map, feasibility framing, figure manifest).

Today's session is a return to the model: tidy the data, recover the original
baseline run, and produce the two mechanistic-finding figures the thesis
Results section claims.

## 1. `out/` tidy-up

Renamed all 30 run directories with a `YYYY-MM-DD_` prefix (date of the files
within) so old data is obvious at a glance. `out/` is gitignored, so no repo
churn; a reference scan confirmed nothing functional pointed at the old names
(only generic example commands, a docstring, historical lab books, and
`reports/figures/v0.5/*.png` filenames).

Identified the thesis-relevant runs and archived the rest (25 exploratory /
superseded runs → `out/_archive/`). Kept at top level:

| Run | Thesis role |
|---|---|
| `2026-05-22_phase3_issue44_final` | most recent validation (v0.5) — Fig. 2 |
| `2026-05-20_ca-bound-free` | current validation-figure source |
| `2026-05-22_dose_sweep_9x9_focus` / `_transition` | dose-response — Fig. 3 |
| `2026-05-12_long-recovery-v041` | DTS overshoot |
| `2026-04-29_original-decay-stub` | regenerated baseline (see §2) |

## 2. Regenerating the original single-protein-decay run

The "model was just a single protein decaying" baseline was never in `out/`
(all runs there are calcium-model runs; `out/` is gitignored, so the original
data was never committed). Recovered it from git history instead:

- `resting_decay.py` was added 2026-04-21 (`768db055`); `calcium_dynamics.py`
  2026-04-29 (`60048186`). The decay-only era is `60048186^`.
- Checked that commit out in an **isolated git worktree** (so the current
  working tree — with all the thesis docs — was untouched), confirmed the April
  code still runs under the current environment, and re-ran
  `verify_resting_decay.py --days 1`.
- Result: proteins decay ~9.4 %/day (half-life 7 days, p = 1.15×10⁻⁶ per 1-s
  step), metabolites unchanged — matches the theoretical 9.43 %. Worktree torn
  down afterwards.

Output: `out/2026-04-29_original-decay-stub/` with a `PROVENANCE.txt` recording
it as a faithful 2026-06-10 re-run, not the literal original.

## 3. PMCA / MCU perturbation experiments (#53)

Built `runscripts/manual/runPerturbation.py` on the `runDoseSweep.py`
service-layer pattern (override one constant → run → harvest → **restore in a
`finally`**), plus a smoke test (`test_perturbation.py`, 2 pass, mypy clean).

### Exploratory run reshaped the design

A quick exploratory run resolved the three open design questions and overturned
the original protocol assumption:

- **Resting state holds** under every perturbation, including `V_max_MCU = 0`
  (cyt 99.9 → 104.9 nM either way) — MCU knockout is safe.
- **The default sustained agonist never returns to rest** — thrombin
  irreversibly cleaves PAR1/4, so cyt rises to a sustained ~436 nM plateau and
  the DTS depletes fully. "Time to recovery" and "DTS refill time" are
  therefore undefined under that protocol. Design refined:
  - **Exp A (PMCA)** → run in the **EDTA / no-extracellular-Ca²⁺** condition,
    where the transient is self-limiting; observable = **recovery-tail AUC**
    (censoring-free).
  - **Exp B (MCU)** → reframed from "DTS refill time" to **cytosolic buffering**.

### Results

- **PMCA Vmax rate-limits the recovery tail.** EDTA transients all peak ~317 nM;
  recovery-tail AUC falls monotonically with Vmax: 46722 → 45668 → 44355 →
  41035 → 35776 nM·s across ×0.25 → ×4. The pump shapes the tail, not the spike.
- **MCU buffers cytosolic Ca²⁺ without rescuing the store.** Peak cyt 651 (KO) →
  436 (baseline) → 380 nM (×4); DTS depletes to 0 in every case. KO reproduces
  the elevated-cyt phenotype of MCU⁻/⁻ platelets (Ghatge 2026). Bonus: ×4 lowers
  the peak but raises the recovery-AUC — slow NCLX release prolongs the tail.

Output: `out/2026-06-10_perturbation/` (`pmca_traces.png`, `mcu_traces.png`,
two `.npz`, `perturbation_summary.json`). The design doc
(`reports/design/perturbation-figures.qmd`) and issue **#53** record the
refinement and tick off the acceptance criteria.

## 4. Wired into the thesis

Copied the two figures (and the decay-stub figure) into `reports/figures/v0.5/`,
embedded them as `@fig-pmca` / `@fig-mcu` in the `draft-thesis.qmd` Results
section, and rewrote the two mechanistic-finding paragraphs from planning notes
into real prose citing the figures. With the validation and dose-sweep figures
this completes the handbook's four-figure Results requirement.

## Files added / changed

- `runscripts/manual/runPerturbation.py` (new, ~250 LoC)
- `models/platelet/tests/sim/test_perturbation.py` (new)
- `reports/design/perturbation-figures.qmd` (refinement note added)
- `reports/figures/v0.5/perturbation-{pmca,mcu}-2026-06-10.png`,
  `original-decay-stub-2026-04-29.png` (new)
- `reports/thesis/draft-thesis.qmd` (Results: 2 figures + prose;
  also corrected a stale `raw_data/` reference in Methods)
- `out/` — 30 dirs renamed, 25 archived, decay-stub + perturbation runs added

## Next steps

- Consider regenerating the validation figure (Fig. 2) from
  `2026-05-22_phase3_issue44_final` for v0.5 consistency (currently cites a
  v0.4.0 figure).
- Insert the dose-sweep figure into Results as Fig. 3 (data already in hand).
- Catalogue the remaining references in Zotero (#52) — incl. Ghatge 2026, now
  cited in Results.
