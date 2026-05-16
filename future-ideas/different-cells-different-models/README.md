# Different cells, different models

A Perspective paper on choosing a modelling paradigm by cell biology, with the human platelet as the worked example.

**Status:** first-draft scaffold (2026-05-16). All sections drafted in prose; quantitative claims and reference fields marked `TODO` where unverified.

**Format:** Quarto + PNAS class. Default target PNAS; may retarget to Cell Systems / npj Systems Biology / PLOS Computational Biology.

## Build

```bash
quarto render paper.qmd --to pnas-pdf
```

Output: `paper.pdf`.

## Layout

| File | Purpose |
|---|---|
| `paper.qmd` | Manuscript source — abstract, significance, body, figure captions |
| `references.bib` | Bibliography (BibLaTeX) — every entry stubbed with `TODO` where fields are unverified |
| `_extensions/christopherkenny/pnas/` | PNAS Quarto template (Christopher Kenny, MIT) |
| `_quarto.yml` | Project config |
| `figures/` | Figure source (TikZ) — to be populated as figures are finalised |

## Authors

- Steve Haigh — School of Biological Sciences, University of Reading
- Dr Mike Fry — University of Reading (placeholder, to be confirmed)

## Sibling work

This Perspective draws on the platelet whole-cell modelling work in
the parent repository (`platelet-wcm`). Figure 3 (platelet case study)
maps directly onto current and planned subsystems of that model.
