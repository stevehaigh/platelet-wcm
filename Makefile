.PHONY: help run stop deploy pdfs pdfs-clean

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
	@echo "    pdfs          Build PDFs from reports/*.md into reports/pdf/"
	@echo "    pdfs-clean    Delete reports/pdf/"
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

