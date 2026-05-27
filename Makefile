.PHONY: help run stop deploy pdfs pdfs-clean kinetics-review kinetics-fonts

help:
	@echo ""
	@echo "  platelet-wcm — available targets"
	@echo ""
	@echo "  Webapp"
	@echo "    run           Run webapp with hot-reload (port $(PORT))"
	@echo "    stop          Kill a running webapp process"
	@echo ""
	@echo "  Azure"
	@echo "    deploy        Push current branch to webapp → triggers Azure CI pipeline"
	@echo ""
	@echo "  Reports"
	@echo "    pdfs            Build PDFs from reports/*.md into reports/pdf/"
	@echo "    pdfs-clean      Delete reports/pdf/"
	@echo "    kinetics-review Render reports/design/kinetics-v0.5-review.pdf from calcium-v0.5.toml"
	@echo ""
	@echo "  Options"
	@echo "    PORT=NNNN     Override port (default: $(PORT))"
	@echo ""

PORT ?= 8050

.DEFAULT_GOAL := help

# ── Local development ─────────────────────────────────────────────────────────

run:
	PYTHONPATH="$$PWD" python3 runscripts/manual/webapp.py --port $(PORT) --debug

stop:
	pkill -f "webapp.py" && echo "Stopped." || echo "No webapp process found."

# ── Azure deployment ──────────────────────────────────────────────────────────

deploy:
	@echo "Pushing $$(git branch --show-current) → webapp to trigger Azure deployment..."
	git push origin HEAD:webapp

# ── Reports: Markdown → PDF ───────────────────────────────────────────────────

REPORTS_MD  := $(shell find reports -name '*.md' -not -path 'reports/pdf/*')
REPORTS_PDF := $(patsubst reports/%.md,reports/pdf/%.pdf,$(REPORTS_MD))

pdfs: $(REPORTS_PDF)
	@echo "Built $(words $(REPORTS_PDF)) PDF(s) into reports/pdf/"

# External-audience reports: no ToC, title comes from the document's YAML
# frontmatter rather than the filename. Listed before the generic rule so
# Make picks the more specific pattern for `reports/external/*.md`.
reports/pdf/external/%.pdf: reports/external/%.md reports/pandoc-header.tex
	@mkdir -p $(dir $@)
	pandoc "$<" -o "$@" \
		--pdf-engine=xelatex \
		--variable=geometry:margin=2.5cm \
		--variable=fontsize:11pt \
		--variable=mainfont:"STIX Two Text" \
		--variable=monofont:"Menlo" \
		--variable=colorlinks:true \
		--variable=linkcolor:blue \
		--variable=urlcolor:blue \
		-H reports/pandoc-header.tex

reports/pdf/%.pdf: reports/%.md reports/pandoc-header.tex
	@mkdir -p $(dir $@)
	pandoc "$<" -o "$@" \
		--pdf-engine=xelatex \
		--variable=geometry:margin=2.5cm \
		--variable=fontsize:11pt \
		--variable=mainfont:"STIX Two Text" \
		--variable=monofont:"Menlo" \
		--variable=colorlinks:true \
		--variable=linkcolor:blue \
		--variable=urlcolor:blue \
		-H reports/pandoc-header.tex \
		--toc --toc-depth=2 \
		--metadata title="$(notdir $(basename $<))"

pdfs-clean:
	rm -rf reports/pdf

# ── Reports: Quarto (.qmd) → PDF ──────────────────────────────────────────────
# Quarto handles native mermaid → PDF rendering without the separate mmdc step.
# Use .qmd for diagram-heavy or interactive-output design docs; plain prose
# stays as .md and goes through the pandoc rule above.
#
# Output goes to reports/pdf-quarto/ so it sits alongside (not on top of) the
# pandoc-built PDFs in reports/pdf/ — gives you a side-by-side comparison.

REPORTS_QMD     := $(shell find reports -name '*.qmd')
REPORTS_QMD_PDF := $(patsubst reports/%.qmd,reports/pdf-quarto/%.pdf,$(REPORTS_QMD))

quarto-pdfs: $(REPORTS_QMD_PDF)
	@echo "Built $(words $(REPORTS_QMD_PDF)) Quarto PDF(s) into reports/pdf-quarto/"

reports/pdf-quarto/%.pdf: reports/%.qmd reports/pandoc-header.tex
	@mkdir -p $(dir $@)
	quarto render "$<" --output-dir "$(CURDIR)/$(dir $@)"

quarto-pdfs-clean:
	rm -rf reports/pdf-quarto

# ── Kinetics review: TOML → clickable PDF ─────────────────────────────────────
# Regenerates reports/design/kinetics-v0.5-review.{qmd,pdf} and
# reports/params/calcium-v0.5-references.bib from
# reports/params/calcium-v0.5.toml. The runscript shells out to `quarto render`
# for the PDF step, so quarto + xelatex must be on PATH.
#
# The QMD header sets `mainfont: "TeX Gyre Termes"` + `monofont: "DejaVu Sans
# Mono"` — both ship with TeX Live but xelatex/fontspec resolves them through
# the OS font lookup (CoreText on macOS, fontconfig on Linux), not TeX's own
# kpathsea. `kinetics-fonts` copies them out of the TeX Live tree into the
# user font dir so fontspec can find them. Idempotent.

KINETICS_FONT_DIR := $(if $(filter Darwin,$(shell uname)),$(HOME)/Library/Fonts,$(HOME)/.fonts)

kinetics-fonts:
	@set -e; \
	 if command -v fc-list >/dev/null 2>&1 \
			&& fc-list | grep -q "TeX Gyre Termes" \
			&& fc-list | grep -q "DejaVu Sans Mono"; then \
		echo "kinetics-fonts: TeX Gyre Termes + DejaVu Sans Mono already installed"; \
		exit 0; \
	 fi; \
	 mkdir -p "$(KINETICS_FONT_DIR)"; \
	 gyre_dir=$$(dirname $$(kpsewhich texgyretermes-regular.otf 2>/dev/null)); \
	 if [ -z "$$gyre_dir" ] || [ ! -d "$$gyre_dir" ]; then \
		echo "error: TeX Gyre Termes OTFs not found via kpsewhich; is TeX Live installed?" >&2; \
		exit 1; \
	 fi; \
	 for f in "$$gyre_dir"/texgyretermes-*.otf; do \
		cp -n "$$f" "$(KINETICS_FONT_DIR)/" 2>/dev/null || true; \
	 done; \
	 dejavu_path=$$(kpsewhich DejaVuSansMono.ttf 2>/dev/null); \
	 if [ -z "$$dejavu_path" ] || [ ! -f "$$dejavu_path" ]; then \
		echo "error: DejaVu Sans Mono TTFs not found via kpsewhich" >&2; \
		exit 1; \
	 fi; \
	 dejavu_dir=$$(dirname "$$dejavu_path"); \
	 for f in "$$dejavu_dir"/DejaVuSansMono*.ttf; do \
		cp -n "$$f" "$(KINETICS_FONT_DIR)/" 2>/dev/null || true; \
	 done; \
	 if command -v fc-cache >/dev/null 2>&1; then fc-cache -f "$(KINETICS_FONT_DIR)" >/dev/null; fi; \
	 echo "kinetics-fonts: TeX Gyre Termes + DejaVu Sans Mono installed under $(KINETICS_FONT_DIR)"

kinetics-review: kinetics-fonts
	PYTHONPATH="$$PWD" python runscripts/manual/buildKineticsReview.py

