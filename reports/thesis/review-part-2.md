# Academic Advisor Review — Platelet WCM MRes Thesis

**Part 2: strategy, thesis-draft critique, and feasibility framing**
Date: 2026-06-09 · Advisory notes for Steve Haigh

---

## 0. Decisions locked this session

These four answers shape every recommendation below.

| # | Question | Decision |
|---|----------|----------|
| 1 | Main-thesis word limit | Guide, not hard. **Target 4,500–5,000 words.** The handbook sub-limits roughly fit this, so no aggressive cutting to supplements is needed. |
| 2 | Supplements assessed? | **No** — optional, almost certainly excluded from the rubric. Wanted anyway: record the work, show engagement/credibility, and double as the PhD plan. Can be brief; **AI may draft substantive parts.** |
| 3 | What does the rubric reward? | **Lab-based focus: methods rigour + solid data analysis.** Feasibility here means *"show the approach is plausible and reason about the gap."* |
| 4 | Claim strength | **Modest.** Show it is *worth pushing the idea further*. Do **not** claim a WCM has been built, or even that it is certainly feasible. Drop any "first mammalian WCM" primacy flag. |

**The through-line:** assessment splits **70 (thesis) / 20 (reflective document) / 10 (presentation)**. The thesis rubric rewards methods and data analysis; the claim is modest; the **reflective document is a second *assessed* home** for the AI-use and decision-making material; the genuinely unassessed supplements are cheap (AI-assisted). Effort priority: **main thesis (esp. Methods + Results) → reflective document → presentation → supplements.**

> **Marking guidance now in hand** (Appendix 4, DISTINCTION descriptors) — mapped to concrete actions in **§9**. It strongly reinforces the Methods + Results emphasis and reveals that **Critical analysis / scientific reasoning is a full, distinction-defining category** in its own right. Graphical-abstract guidance (Jambor et al. + the five concept drafts) is in **§8**.

---

## 1. Where the project stands

- **Biology is frozen at v0.4.1** (2026-05-12). That was the last functional change. The "v0.5" label is the kinetics-as-data refactor (#32: ~150 constants → `calcium-v0.5.toml`, 63 species → `species-v0.5.tsv`), the 2-D dose sweep (#45), and its webapp tab (#46) — **infrastructure and reportability, not new mechanisms.**
- **Healthy:** 52/52 tests pass, mypy clean.
- **The science:** 12 mechanisms, agonist → Ca²⁺ end-to-end, no hand-fitted forcing, Dolan & Diamond 2014 Fig. 4 validated 5/5, resting state held (cyt 104 nM, DTS 235 µM, IP₃ 50 nM).
- **Gap in the written record:** the lab-book trail goes cold at 2026-05-22. The #32 refactor, the TOML reference-verification pass, and the entire thesis draft (June 1–9) have no lab-book entry. Worth a short catch-up note so the work is traceable.

---

## 2. Headline recommendation — write, don't build

**Stop net-new biology. Write.** With one coding task and one strategic caveat.

This is a **feasibility study**, and decision (3) confirms the bar is *"show the approach is plausible and reason about the gap."* You have already done enough simulation to clear that bar: one validated pathway in an extensible, mass-balanced framework. Adding a third receptor or a granule process would strengthen the "extensible" claim only marginally and would not change the feasibility verdict — while trading directly against the writing that earns the degree. Writeup is late June; today is the 9th. Time is the binding constraint.

### The one coding task worth doing (rubric-aligned)

**Generate the perturbation figures your Results section already claims.** The draft asserts "PMCA Vmax rate-limits recovery" and "MCU buffering does not accelerate DTS refill" with **no figure or data behind either**. Because the rubric rewards *data analysis*, these are not optional polish — they are exactly the kind of result you are marked on. Both are single-parameter scans / a knockout run on the frozen v0.5 model (~2–4 days, low risk). Producing **both** gives two genuine data figures and comfortably fills the four-figure Results quota (see §6).

### The strategic caveat (decided: don't build it now)

The one thing that would genuinely *upgrade* the feasibility proof is a **minimal second process sharing the mass-balanced state** (e.g. a granule-release stub that consumes cytosolic Ca²⁺). Today, the "extensible" claim is **asserted, not demonstrated** — you have one pathway, not two composing over one shared state, which is the whole reason you chose wcEcoli over COPASI/VCell. Given decision (4) (modest claim) and the time available, **do not half-build this.** Name it in the Discussion as *the decisive next experiment* — which is exactly the right move for a "worth pushing further" thesis.

---

## 2A. Assessment split (70 / 20 / 10) — and how it redirects effort

The thesis is **70 marks**, but it is not the whole assessment:

| Component | Marks | Implication |
|---|---:|---|
| **Thesis** (~4,500–5,000 words) | **70** | Primary effort. Methods + Results are rubric-weighted (§6). |
| **Reflective document** | **20** | A *second assessed home* — worth 2× the presentation. Much of the "supplement" material (AI use, the k₃ episode, the framework-choice decision, what the work taught you) belongs **here**, where it earns marks, not in unassessed supplements. |
| **Presentation** | **10** | Smallest component; cheap once the thesis figures are locked — reuse them. |

**The key consequence: the reflective document changes the supplements calculus.** §5 below treated the AI-use and decision-making narrative as unassessed supplement material. It is not — it is exactly what a reflective document rewards. Redirect it:

- **Into the reflective document (assessed, 20 marks):** how you used AI and where it failed (the Purvis k₃ reciprocal-error story is a *perfect* reflective episode — AI confident, wrong, caught by primary-source discipline); the framework-choice decision (wcEcoli vs COPASI/VCell) as a reasoned judgement under uncertainty; what building one pathway taught you about the data bottleneck; how your working practice changed.
- **Stays as thesis / supplement:** the CZI/Allen *field positioning* (Discussion / lit-review material, not personal reflection); the SWE-for-biologists model description; the detailed module-by-module work breakdown.

So the "AI + CZI/Allen" idea naturally **splits**: the *state-of-the-field / two-paradigms* half feeds the thesis Discussion (or a lit supplement); the *how-I-personally-used-AI* half is reflective-document gold.

**Reflective-document craft.** Use a named framework — Gibbs' reflective cycle, Kolb's experiential-learning cycle, or Schön's reflection-in/on-action — markers reward a recognised structure over free-form musing. Anchor each reflection to a concrete episode (k₃; the dose-sweep "binary store-release" surprise that overturned your expectation; the argmax rise-time artefact you had to diagnose). Be honest about what went wrong — reflective marks come from *insight into the process*, not from a polished success story.

**Presentation (10 marks).** Lowest marginal cost. Backbone = the four thesis figures + the tiered feasibility verdict (§4). Build it *after* the figures are locked; near-zero extra asset creation.

---

## 3. The thesis draft — advisor critique

### Strengths (keep, don't touch)

- **Feasibility framing is correct** and the abstract states it cleanly.
- **The Purvis k₃ story is the standout.** It demonstrates genuine rigour (cross-checking a review against its primary source, finding a reciprocal transcription error) and doubles as honest evidence for the AI-assisted-extraction discussion. Foreground it in Methods. Given the rubric's methods focus, this is one of your highest-value paragraphs.
- **The framework-choice comparison and the COPASI/VCell supplement** are mature, dispassionate, and well-cited. The "field's own verdict" device (Goldberg/Karr) is exactly how you defend a tool choice to an examiner.

### Concerns (priority order)

**1. Soften the claim (decision 4).** The abstract currently ends *"…with further work a whole-cell model is indeed feasible."* That is stronger than you want. Suggested register:

> *"…demonstrated that a major pathway can be modelled accurately within a whole-cell framework, and that the approach is extensible in principle. We do not claim to have shown that a complete whole-cell platelet model is feasible; rather, that the proof-of-concept is encouraging, and the principal obstacle — data availability, not computation — sufficiently well-characterised, to justify pushing the idea further."*

**2. Drop the primacy flag (decision 4).** The line *"no published whole-cell models of any mammalian cells"* is defensible only with a precise definition (*mechanistic / dynamical / mass-balanced*, the Karr→Covert lineage), because an examiner can point at the **Allen Integrated Cell** (a statistical/ML model of cell morphology) as an apparent counterexample. Since you do not want to over-egg it, **demote this from a claim to context**: note that no published *mechanistic* whole-cell model of a mammalian cell yet exists, cite the Covert lab H1 hESC effort as "first step toward" one, and move on. No flag to defend.

**3. Four-figure Results gap.** Currently: validation (1), dose-sweep (1, *not yet inserted into the draft*), long-recovery/DTS-overshoot (1, figure exists). The architecture diagram is a *schematic*, not a data result — with a data-focused rubric, lean on genuine data figures. The perturbation task (§2) closes the gap.

**4. Insert the dose-sweep — it is your strongest novel result and is currently absent from Results.** The "graded upstream (IP₃ spans 50→305 nM over five orders of magnitude) → binary at store release (peak Ca²⁺ flat ~436 nM) → graded again in the integrated response (Ca²⁺ AUC recovers 8× range)" story is precisely the kind of insight that only emerges from an all-components-present model. It is also clean data analysis — directly rubric-aligned.

**5. Figure version mismatch.** ✅ RESOLVED 2026-06-10 — draft now cites `ca-bound-free-edta-2026-06-10.png` (EDTA transient, fixed endogenous-IP₃ model); `v0.4.0` and `2026-05-20` are superseded.

---

## 4. How strongly is feasibility "proved"? — the intellectual core

You asked for help framing this. Cast the proof in **three tiers**; this is both an honest self-assessment and a ready-made Discussion structure. It maps directly onto decision (3): *show plausible, reason about the gap.*

| Tier | Claim | Status |
|------|-------|--------|
| **Weak** | "A platelet Ca²⁺ pathway can be modelled mechanistically in a WCM framework and validated against experiment." | **Solidly proved.** Dolan 5/5, resting state held, no hand-fitted forcing. |
| **Medium** | "The wcEcoli architecture ports to mammalian/platelet biology and is extensible." | **Partially shown.** ✓ for *adaptation* (pruned + rebuilt); ✗ for *extensibility*, which is **asserted, not demonstrated** — one pathway, not two composing over shared state. This is the honest soft spot. |
| **Strong** | "A genuine whole-cell platelet model is feasible." | **Not proved — and per decision (4), not claimed.** Defensible instead: the approach is a viable proof-of-concept, and the principal risk is *data availability, not architecture or computation* — which the field's own blueprint (Szigeti et al.) independently confirms. |

### Rough work estimate to a "genuine" platelet WCM

Order-of-magnitude, with the dominant cost flagged explicitly:

- **Tractability advantage: the platelet is anucleate.** No transcription, translation, replication, or division — exactly where the Karr 2012 *Mycoplasma* WCM (28 submodels, multi-year team) spent the most effort. This is your strongest feasibility argument.
- **Still required (minimally):** the Ca²⁺ signalling core (✓), metabolism/mitochondrial ATP + glycolysis, cytoskeletal/shape-change dynamics, granule biogenesis + secretion, integrin αIIbβ3 inside-out/outside-in, full GPCR receptor complement, membrane lipid/PI dynamics, procoagulant (PS-exposure) machinery — on the order of **8–15 major process modules**.
- **The bottleneck is data, not code.** You built one module in ~4–6 months part-time; the limiting cost was curation/calibration (the k₃ episode is the proof). AI-assisted coding compresses the *engineering*, not the *data* cost.
- **Defensible estimate:** a small team (2–4 people) over **~3–5 years**, with the explicit qualifier that **2–3 subsystems may be currently infeasible without new wet-lab measurement.** Honest verdict: *"feasible in principle for most modules; partially blocked on data today."* Present as person-years dominated by data, not as a coding schedule.

This tiered framing *is* your answer to the rubric: it shows the approach is plausible (weak tier proved) and reasons carefully about the gap (medium/strong tiers + work estimate).

---

## 5. The supplements — reframed for decision (2)

First, per §2A: the *AI-use and decision-making* material has an **assessed** home (the reflective document, 20 marks) — route it there, not here. What remains below is the genuinely unassessed remainder. Because those supplements are **optional and AI-draftable**, the calculus is no longer "which to cut." It is: **produce them cheaply (AI-assisted, human-verified), because they cost little, signal engagement, and form your PhD plan.** Rank by *human attention required*, not by whether to write:

1. **Feasibility evaluation + work estimate** — most human judgment, most PhD-relevant. The *verdict* belongs in the main Discussion (don't bury the punchline in an appendix); the detailed module-by-module breakdown goes here.
2. **Model description (software-engineer's view, for biologists)** — the artefact that records the actual work. Watch the register: motivate every technical choice *biologically*, or it reads as code docs. AI can draft from the codebase + CLAUDE.md; you verify.
3. **AI + CZI/Allen** — distinctive and current. Frame it as a genuine **two-paradigm contrast**, not a survey:
   - *Mechanistic / bottom-up* (Karr → Covert wcEcoli → **your platelet**): equations, parameters, interpretable, causally interrogable — but **data-bottlenecked**.
   - *AI / data-driven / top-down*: CZI's Virtual Cells Platform (NVIDIA partnership, Oct 2025; petabyte-scale), the rBio reasoning model (Aug 2025), the Allen Integrated Cell — predictive and vastly better-resourced, but **correlative, not mechanistic**.
   - The argument: which reaches a "genuine" WCM first, and are they complementary? Your own first-hand AI use (the k₃ story) is a unique data point bridging the two — connect it.
4. **Literature review** — most safely AI-assisted (with citation verification). Ensure it extends, not duplicates, the Intro's mini-review.

> **⚠ Caution if AI drafts supplements.** Your thesis literally contains the cautionary tale (Purvis k₃: AI confidently reproduced a wrong value from a review). Apply the same discipline: **every citation and every numeric claim AI produces must be verified against the primary source.** And disclose supplement AI use honestly in the thesis AI statement — it is on-theme, not a liability.

---

## 6. Concrete writing plan

### Word budget (target 4,500–5,000)

| Section | Target | Notes |
|---------|-------:|-------|
| Abstract | ~250 | Structured; apply the softened claim (§3.1). |
| Significance | ~120 | Non-specialist. |
| Introduction | ~600 | Mini lit review + gap + the *modest* aim. |
| **Methods** | **~800** | **Rubric-weighted.** Framework adaptation; per-mechanism methodology; the k₃ story; reproducibility (tests/CI/Docker). |
| **Results** | **~1,500** | **Rubric-weighted.** Four figures + analysis (see manifest). Describe without interpretation. |
| Discussion | ~1,600 | Feasibility verdict (tiered, §4); limitations; data-bottleneck argument; work-estimate summary; "worth pushing further" close. |
| Conclusion | ~150 | Modest register. |
| **Total** | **~5,020** | Trim Discussion ~50–100 to land inside 4,500–5,000. |

### Figure manifest

| Fig | Content | Asset | Status |
|-----|---------|-------|--------|
| 1 | Architecture schematic | `model-status-graphviz.png` | ✅ exists (schematic — not a "result"; BioRender redraw optional) |
| 2 | Validation vs Dolan Fig. 4 | `ca-bound-free-edta-2026-06-10.png` + phase3 numbers | ✅ done — draft ref updated 2026-06-10 |
| 3 | Dose-response (graded→binary→graded) | `dose-sweep-9x9-panel.png` | ✅ exists; **not yet in draft — insert** |
| 4 | Perturbation: PMCA Vmax scan → recovery time | — | ❌ **needs run** (the one coding task) |
| 4b | MCU knockout → DTS refill | — | ❌ **needs run**; gives a clean 4th data figure |

Because Fig 1 is a schematic and the rubric is data-focused, generating **both** 4 and 4b yields three solid data figures (2, 3, 4) plus a fourth (4b) — comfortably satisfying the requirement with analysis you are actually marked on.

### Suggested sequence

1. **Lock the figures first** (≤1 week): refresh the validation ref, insert the dose-sweep, run the PMCA/MCU perturbations. Results cannot be written before the figures exist.
2. **Draft Methods + Results** (the marked, human-effort core).
3. **Draft Discussion** around the tiered feasibility framing (§4).
4. **AI-assist the supplements** in parallel/after, human-verifying every citation and number.
5. **Housekeeping:** catch-up lab book (§1), finish the Zotero bibliography (#52).

---

## 7. Quick cautions

- **Lab-book gap** since 2026-05-22 — write the catch-up note while the #32/thesis work is fresh.
- **Bibliography (#52)** is still open and blocks a clean reference list.
- **AI citation discipline** — see §5 warning; your own thesis is the case study.

---

## 8. Graphical / visual abstract

**Source caveat.** The referred paper (Jambor et al. 2021, *PLOS Biology*, "Creating clear and informative image-based figures") is about image-based figures generally — and is itself wet-lab oriented (scale bars, fluorescence channels, micrographs). Take its *transferable* principles and pair them with standard graphical-abstract rules.

**Principles that transfer (Jambor + graphical-abstract convention):**

- **One take-home message.** Minimal text; self-explanatory to a non-specialist (the significance-statement audience).
- **Clear reading direction** (top→bottom or left→right) with white space separating blocks.
- **Colourblind-safe, perceptually-uniform colormap** (viridis or magma) used *consistently across all panels*; never convey meaning by colour alone — keep axis labels and values. (Jambor: ~45% of cell-biology papers had deuteranopia-inaccessible images.)
- **Units + range on every quantitative element.** Jambor's "every image needs a scale bar" → for you, **every heatmap needs a colourbar with units (nM, nM·s) and a numeric range.** The current wireframes have none — this is the single biggest gap.
- **Define or drop abbreviations** (DTS, SOCE, AUC) on the figure itself.

**The five concept drafts:**

| Concept | Message | Strength | Risk | Best role |
|---|---|---|---|---|
| 1 — Three-layer signal shaping | graded → binary → graded | Encodes the dose-sweep story directly | 6 boxes + 3 heatmaps = dense for a GA | Dose-sweep **Results figure** |
| 2 — Single-cell biology constrains dose-response | platelet + one inset | Compact | Subtle message; overlapping arrow labels; needs domain knowledge | Weakest as a GA |
| 3 — WCM funnel (12 mechanisms) | comprehensive platform | Shows scope | Jargon-list-heavy; diffuse, not a single insight | Architecture/scope figure |
| 4 — Why peak Ca²⁺ hides the graded response | peak binary vs AUC graded | **Sharpest single insight; most accessible**; drug-screening hook | Carries a *finding*, not the feasibility *contribution* | Strong GA **or** headline Results figure |
| 5 — Prior single-pathway → whole-cell platform | from one ODE to integrated cell | **Visualises the feasibility contribution** (the thesis's actual punchline) | Jargon in layer labels; needs simplifying for a broad audience | Strong GA for a *feasibility* thesis |

**Recommendation.** A feasibility thesis's message is *"the whole-cell approach is worth pursuing,"* so the graphical abstract should sell the **contribution** → **Concept 5 is the best fit**, but dial the jargon back (the *number* and *integration* of mechanisms matter more than `MCU/NCLX/CaM/buffers`; use a simplified stacked cartoon). **Strong alternative: Concept 4** if you would rather lead with the headline finding — it is the most accessible and striking of the five. Then reuse the rest as body figures (Concept 1 or 4 for the dose-sweep result; Concept 3/5 structure for the architecture figure).

**Universal fixes for whichever you pick:** add colourbars with units + range to every heatmap; standardise on one colourblind-safe colormap (you are currently mixing a viridis-like and a magma-like map); replace wireframe boxes with the planned simplified BioRender cartoons; keep total text under ~30 words; render at ≥300 dpi with legible fonts.

---

## 9. Marking descriptors (DISTINCTION) → concrete actions

**Your caveat is correct and load-bearing:** these descriptors are written for *wet-lab* projects. Two bullets ("appropriate controls" and "use of statistics") do not map cleanly onto a deterministic in-silico model. **The distinction-level response is not to ignore them but to translate them explicitly and justify the translation** — which itself satisfies the Critical-Analysis bullet *"reasoned reflection of methodology and questioning of assumptions."* Engaging with the mismatch earns marks; silently skipping it loses them.

| Descriptor category | What it rewards | Your action / asset |
|---|---|---|
| **Presentation** | Layout, clarity, correct terminology + **units throughout**, spelling/grammar, **consistently formatted references** | Quarto + Harvard CSL already in place; **finish the bibliography (#52)** — "appropriate, consistently formatted references" is an explicit bullet. Enforce unit discipline (µM, nM, s). |
| **Introduction** | Mastery of **wide literature**, **wider reading / beyond standard material**, clear aims/hypothesis | The CZI/Allen field-positioning earns the **"wider reading / beyond standard material"** mark *here* (not only in the reflective doc). State a crisp aim: a feasibility study's "hypothesis" is *"a major platelet pathway can be modelled mechanistically and validated within a whole-cell framework."* |
| **Methods** | Clearly written, **repeatable**, grasp of techniques, work suitable/**exceeds** the duration | Your reproducibility story is *stronger* than most wet-lab repeatability: version-controlled code, pinned deps, fixed seeds, `OPENBLAS_NUM_THREADS=1`, CI, Docker, public repo. Make it explicit. The 12-mechanism scope + dose-sweep evidences "exceeds the work expected." |
| **Results** | Exemplary figures + **concise legends**; **statistics where applicable**; **controls where applicable**; design appropriate to the hypothesis; significant achievement | See translations below. §6 figure manifest + §8 cover figures/legends. |
| **Critical analysis / scientific reasoning** *(the largest, distinction-defining category — 7 bullets)* | Integration of theory + evidence; balanced argument; **reflection on methodology + questioning assumptions**; lateral thinking/insight; **critique + counter-critique**; independent judgement; **mature limitations + future work** | This is where a distinction is won or lost, and where your honest style shines. Feeds the Discussion **and** the reflective document. Assets: the tiered feasibility framing (§4); honest limitations (SERCA Vmax over-estimate, resting-baseline tension #37/#38, well-mixed / no microdomains, the DyK low-Ca²⁺ bistability artefact); the COPASI/VCell "roads not taken" counter-critique; the k₃ story (= questioning assumptions); the work estimate (= future work). **Do not over-claim** (decision 4) — maturity here means calibrated confidence. |

**The two wet-lab bullets, translated:**

- **"Excellent inclusion and use of appropriate controls (where applicable)."** Name your *in-silico* controls explicitly so the marker sees you engaged with the concept: the **Dolan ±extracellular-Ca²⁺ paired condition** (the validation control), the **resting / at-rest baseline**, and the **knockout / parameter-scan perturbations** (MCU-off, PMCA-Vmax scan) as controlled manipulations. These are genuine controls — label them as such.
- **"Mature and fully justified use of statistics to analyse the data (where applicable)."** A deterministic ODE has no sampling distribution, so classical inferential statistics are largely **not applicable** — *and saying so, with the reasoning, is the mature move the descriptor rewards.* The legitimate analogue is **sensitivity analysis** (vary a parameter, report the effect on observables) and the **dose-sweep parameter scan** you already have. Frame these as your quantitative analysis; state plainly why error bars / p-values are not the right tool for a deterministic model.

---

## 10. Division of labour across the document set

The whole-cell lineage (Karr → Covert → Szigeti) and the "bottleneck is data, not computation" argument recur across several documents. To avoid repetition (which an examiner will notice), give each its own job and cross-reference rather than restate:

| Document | Owns | Defers to |
|---|---|---|
| Main thesis (Intro + Discussion) | the headline argument + the feasibility verdict (brief) | feasibility supplement for the evidence |
| Lit review | the lineage as *literature* (what has been done) | AI supplement for the paradigm contrast |
| AI supplement | the lineage as *one of two paradigms* (mechanistic vs data-driven) | reflective doc for the personal AI account |
| Feasibility supplement | the *evidence + module-by-module work estimate* | main Discussion for the one-line verdict |
| Reflective document | the *personal* process narrative (k₃, decisions, development) | — |

State the data-bottleneck thesis **in full once** (feasibility supplement); reference it elsewhere.

**Word budgets (supplements are unassessed — don't over-invest):** cap lit review, model description, and AI supplements at **~1,500 words each**; the feasibility supplement may run to **~2,500** (it carries the most original reasoning). The reflective document is **~2,000 and assessed** — that is where writing time goes after the main thesis.

---

## Sources

- [Allen Integrated Cell](https://www.allencell.org/allen-integrated-cell.html)
- [A blueprint for human whole-cell modeling — Szigeti et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5966287/)
- [CZI–NVIDIA Virtual Cells Platform (Oct 2025)](https://chanzuckerberg.com/newsroom/nvidia-partnership-virtual-cell-model/)
- [CZI rBio reasoning model (Aug 2025)](https://chanzuckerberg.com/blog/rbio-reasoning-ai-model/)
- [Jambor et al. 2021, *PLOS Biology* — Creating clear and informative image-based figures](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001161)

*Next (if wanted): section-by-section line edits of the draft once the figures are locked, and/or a detailed module breakdown for the feasibility-evaluation supplement.*
