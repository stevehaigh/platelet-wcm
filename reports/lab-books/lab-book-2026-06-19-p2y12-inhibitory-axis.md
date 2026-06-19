---
title: "Lab book — 2026-06-19: the inhibitory axis (P2Y12/Gi/cAMP/PKA), the granule-secretion plot, and a validation re-think"
---

# Lab book — 2026-06-19: P2Y12 inhibitory axis, secretion plot, validation reframe

## Context

A long, productive session that closed two issues and re-framed how the model is
validated. Three threads: (1) a small framework-integration job — the
granule-secretion analysis plot (#8); (2) the substantive science — adding
**P2Y12 / Gi / cAMP / PKA**, the first piece of the v0.7 *inhibitory* axis and
the clopidogrel/ticagrelor drug target (#10); and (3) a deliberate re-think of
the "Dolan 5/5" validation now that the model has outgrown its Ca²⁺ core. Both
issues shipped (PRs #66 and #68 merged to `main`); an `xhigh` multi-agent code
review ran over #68 and its findings were applied before merge. Notes captured
here because the work spanned two PRs, a literature-mining pass, a design
discussion, a figure rebuild, and a 15-point review.

---

## 1. Granule-secretion kinetics plot — #8 (PR #66)

The secretion *data* already existed (v0.61 `SecretionTrace`); what was missing
was a framework-integrated plot, as opposed to the standalone
`plotDownstreamModules.py secretion_release` figure.

- Added `models/platelet/analysis/single/granule_secretion.py` as a
  `SingleAnalysisPlot` reading `SecretionTrace` + `CalciumTrace`: a 3-panel
  figure (cytosolic Ca²⁺ trigger; cargo released / surface-exposed %; cumulative
  secreted autocrine ADP + serotonin).
- Registered in `single/__init__.py` `ACTIVE`, so it runs via
  `analysisPlatelet.py --plot granule_secretion` **and** appears in the webapp
  Explore tab (the webapp runs the default analysis set; `find_plot_images`
  globs the resulting PNG).
- Verified the release curve is sigmoidal (ADP 0 % at rest → 67 % at t=50 s →
  99.9 % plateau; secretion gate 0 → 0.37). Onset uses threshold-crossing, not
  `argmax` (the curves are flat-topped).
- Tests added; the default-plot-list assertion updated. Merged squash `f3f1f8d7`.

A scope note baked into the docstring: granules are bulk cargo pools, so release
is shown as the fraction of each conserved pool that reached `[e]`; the
per-instance `loaded → fusing → released` view still awaits the UniqueMolecule
refactor (#5, deliberately deferred this session).

---

## 2. The inhibitory axis — P2Y12 / Gi / cAMP / PKA — #10 (PR #68)

The substantive science. Through v0.6/v0.61 the model became a rich *activation*
model (Gq cascade → PI cycle → Ca²⁺ core, PKC feedback, secretion / thromboxane
/ integrin outputs) with **no representation of inhibition**. #10 adds ADP's
second receptor — **P2Y12 (Gi-coupled)** — and the cAMP/PKA "off" node it acts
on. Biology: `ADP → P2Y12 → Gi → ↓adenylyl cyclase → ↓cAMP → ↓PKA →
dis-inhibits the cascade`. P2Y12 reads the *same* pericellular ADP as P2Y1 (one
ligand, two receptors) but couples to Gi, not Gq.

### Scope decision

A complete v0.7 inhibitory-axis design already existed
(`reports/design/inhibitory-axis-design-2026-06-15.qmd`), which sequences P2Y12
as **Slice 2** — and crucially makes it *depend on* a cAMP/PKA core (Slice 1):
"ADP lowers cAMP → less PKA brake" is meaningless without a PKA brake to relieve.
So P2Y12 cannot ship alone. Chosen scope (after surfacing the trade-off): a
**lean** cAMP/PKA core + P2Y12 + the VASP/PRI readout, deferring PGI₂/Gs, an
explicit PDE3A, and the NO/cGMP/PKG arm to later slices.

### Parameters (literature-mined, verified)

- **P2Y12: 425 ± 51 sites/platelet, Kd 3.3 nM** — [³H]PSB-0413 radioligand
  binding, reported in Gachet 2012 (*Purinergic Signal.* 8:609; DOI verified via
  the PMC page). Added as `[references.gachet_2012]`. P2Y1 is ~150 (matches the
  existing model).
- **Resting cAMP ~1 µM**; PKA half-max at cAMP ~0.1–0.3 µM; platelets default to
  elevated cAMP / active PKA (the resting quiescence tone).
- **VASP ~45 000–79 000/platelet** (Burkhart 2012); used 50 000.

### Implementation (5 new species + algebraic PKA)

- `P2Y12_inactive/active[pl]`, `cAMP[c]`, `VASP[c]/VASP_phos[c]` in
  `species-v0.6.tsv` + `MOLECULE_NAMES`.
- ODE (`calcium_signalling._ode_rhs`): P2Y12 reversible ADP binding (**not** in
  `total_active_R` — it's Gi); cAMP node
  `d/dt = V_AC·(1 − i_gi_max·p2y12_frac) − k_pde·cAMP`; PKA = Hill(cAMP),
  algebraic; a shared `pka_brake_factor()`.
- **The load-bearing trick**: the PKA brake is `b = 1 + gain·(PKA_rest − PKA)`,
  *normalised to exactly 1.0 at resting cAMP*. Resting ADP = 0 → no P2Y12
  occupancy → cAMP at basal → b = 1, so the resting fixed point and Dolan 5/5 are
  preserved **by construction**, and P2Y12 can only amplify above baseline. (The
  review later made this exact — see §6.)
- `RunConfig.p2y12_block ∈ [0,1]` is the competitive-antagonist knob
  (cangrelor/ticagrelor); registered in the TUI RunSpec schema and added to the
  knockout entity map.
- `CalciumTrace` gained `camp_uM`, `pka_active_frac`, `p2y12_active_frac`,
  `vasp_phos_frac`.

Merged squash `c7afdf65`; closes #10. Dry-mass baseline shifted 261.70 → 265.05
fg (the VASP pool); byte-identical goldens regenerated; full suite green.

---

## 3. The key finding — cytosolic Ca²⁺ is architecturally clamped

The IP₃R PKA brake **does essentially nothing to cytosolic Ca²⁺**. Probes
(while calibrating): IP₃R open-prob brake gain 1→8 moved the peak 293→298 nM;
PMCA k_cat 1.0→0.3 moved 408→412 nM; SERCA 1.0→0.3× moved ~0. The cytosolic
level is pinned by the buffers (CaM/gelsolin/calreticulin) + the SOCE↔clearance
balance, **independent of any single release/clearance rate**. This generalises
the earlier `pkc_ca_invisible` finding from "PKC feedback is invisible" to "*any*
upstream/clearance lever is invisible on free cytosolic Ca²⁺ without a deep
Ca²⁺-core recalibration (which would threaten Dolan 5/5)."

**Resolution — brake the output, not the cytosol.** P2Y12's biologically
dominant effect (and clopidogrel's actual mechanism) isn't cytosolic Ca²⁺ at
all — it's amplifying **integrin activation / secretion / sustained
aggregation**. So the *visible* PKA brake was put on the **αIIbβ3 inside-out
rate** (`IntegrinActivation`), where it is not Ca²⁺-clamped and is biologically
correct. Result: clopidogrel (`p2y12_block=1`) cuts PAC-1 activation **83 % →
70 %** (standard agonist) while keeping cAMP/VASP-P high — the real drug
mechanism, visible. The IP₃R brake is kept too (real, Dolan-safe, but small).

This is the general escape hatch for any "invisible-on-Ca²⁺" feedback in this
model: brake the functional output (integrin/secretion/thromboxane), not the
clamped cytosolic level.

---

## 4. Antiplatelet treatments figure — and the aspirin question

Built `runscripts/manual/plotInhibitoryAxis.py` → two figures in
`reports/figures/v0.7/`:

- `inhibitory_axis_mechanism.png` — ADP → P2Y12 occupancy drives cAMP down from
  ~1 µM; PKA and the VASP/PRI readout fall and recover with it.
- `antiplatelet_treatments.png` — control vs aspirin vs clopidogrel vs both.

**The aspirin puzzle.** First cut (saturating agonist) made aspirin look inert.
Trace confirmed it is *not*: aspirin (`cox1_factor=0`) zeroes TXA₂ and, in a
non-saturated regime, cuts active Gαq **4683 → 3624** by removing the autocrine
TXA₂→TP→Gαq loop. It only *looked* dead because at a saturating agonist Gαq is
already pinned at its ~5000 pool ceiling (so removing the TXA₂ contribution
changes nothing downstream) — clinically reasonable, since aspirin is a weak
antiplatelet against strong thrombin. Rebuilt the figure in the **weak-agonist
(autocrine second-wave) regime** with a 6th panel (active Gαq), so each drug is
visible on its own lever: clopidogrel via cAMP/PKA → integrin (A,C,D); aspirin
via TXA₂ → Gαq (B,E); cytosolic Ca²⁺ (F) is the shared, clamped output.

---

## 5. Validation re-think — "Dolan 5/5" demoted to a regression invariant

Prompted by a good question: is "5/5" still valuable? Conclusion — **no longer
as the headline**. Three reasons, each established this session: (a) it now
passes *by construction* (every new layer is normalised at rest to keep it
passing), so it guards rather than validates; (b) the Ca²⁺ clamp (§3) gives it
low discriminating power — many parameterisations pass it; (c) the model has
outgrown the Ca²⁺ transient (PI cycle, GPCR cascade, PKC, secretion,
thromboxane, integrin, now the inhibitory axis).

Action: wrote `reports/design/validation-map-2026-06-19.qmd` — reframes 5/5 as a
**regression invariant on the Ca²⁺ core**, lays out a *portfolio* of
subsystem-specific targets (Ca²⁺ → Dolan; integrin → PAC-1 flow cytometry;
inhibitory axis → VASP/PRI; secretion → lumi-aggregometry; drug dose-response),
and documents the Ca²⁺-clamp limitation as honest future work. The
`test_phase3.py` 5/5 docstring was reworded to point at this. This session's work
already added two genuine portfolio targets (VASP/PRI + a PAC-1 drug-response)
that the Ca²⁺ transient could never provide.

---

## 6. Code review (`xhigh`, multi-agent) + fixes

Ran the multi-angle review over the #68 diff (10 finder angles → verify → sweep);
15 findings, all applied before merge. Highlights:

- **Resting brake made *exactly* 1.0.** The review caught that the cAMP initial
  count (integer 3613) disagreed with the AC/PDE fixed point (3613.2), so b was
  1.0000084 at t=0 (a tiny self-correcting transient that perturbed the at-rest
  DTS golden). Fixed by `CAMP_REST_COUNT = round(...)` and deriving V_AC + the
  brake anchor from it → `b_ip3r(3613) = b_integrin(3613) = 1.0` exactly. The "by
  construction" claim is now literally true.
- **`p2y12_block` clamped to [0,1]** (block > 1 had flipped the ADP on-rate
  negative → non-physical inverted drug response, no error); AC term floored at 0
  and `p2y12_frac ≤ 1` so a mis-set `i_gi_max > 1` can't drive cAMP production
  negative.
- **Dead `b_min` override kwarg removed**; `pka_active_frac` computed once per ODE
  substep (was twice); `PKA_REST_FRAC` reuses the helper.
- **`plotInhibitoryAxis` reads listener columns** (`IntegrinTrace.active_frac`,
  `ThromboxaneTrace.txa2_uM`) and opens each table once per run (was 3–4 full
  BulkMolecules reads/run); αIIbβ3 → mathtext; unused import dropped; a
  function/attribute name collision renamed.
- **Integrin PKA brake now logged** (`IntegrinTrace.pka_brake`) — the clopidogrel
  mechanism was previously recorded nowhere; **P2Y12 added to
  `knockouts.ENTITIES`** (appended, to keep the TUI's index-based checkbox tests
  stable) for a uniform expression-knockout story.

Final: platelet **122** + wholecell **117** tests pass, mypy clean, CI green
across pytest / mypy / kinetics-review.

---

## 7. Housekeeping & notes

- **Branch cleanup**: merged branches removed; `main` synced. (Two non-deploy
  branches — `model-diagram-checklist`, `claude/platelet-web-agent-model-…` —
  left for manual handling per earlier instruction; `webapp` is the deploy
  branch, kept.)
- **Copilot code-review keeps failing** — diagnosed: it is pinned (GitHub-side) to
  model `claude-opus-4.8`, which the Copilot review backend reports as *not
  available*, so `session.create` fails before reviewing. Per GitHub docs the
  code-review model is **not user-selectable** ("Model switching is not
  supported"; the Models settings page governs Chat only). It is non-blocking and
  unrelated to the code; our real CI gates the PR. Nothing to change in-repo.

## State / next

- v0.7 **Slice 2 (P2Y12/Gi/cAMP/PKA) shipped**; #8 and #10 closed.
- Remaining inhibitory-axis slices (documented, not yet built): full Slice 1
  (PGI₂→Gs→AC + explicit PDE3A + cAMP-elevating drug knobs: iloprost/forskolin/
  cilostazol); Slice 3 (NO/cGMP/PKG); Slice 4 (PKA brake on the secretion-side
  output). A re-calibration toward a *mechanistic* resting cAMP tone (vs the
  current perturbative, normalised-at-rest brake) is the deeper future step — and
  the natural place to revisit the Ca²⁺-core buffer/SOCE clamp from §3.
