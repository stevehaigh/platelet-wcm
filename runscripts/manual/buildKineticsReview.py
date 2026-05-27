"""buildKineticsReview.py — render a clickable PDF review document
from `reports/params/calcium-v0.5.toml`.

Walks the TOML file line-by-line to recover the comment structure
(which `tomllib` discards), groups it by section, and emits a Quarto
``.qmd`` source. Citations of the form ``Author YEAR`` are auto-linked
to entries in the TOML's ``[references.*]`` section; bare DOIs (``doi:X``)
and URLs (``https://...``) are also linkified.

Usage
-----
    PYTHONPATH=$PWD python runscripts/manual/buildKineticsReview.py

Output
------
- `reports/design/kinetics-v0.5-review.qmd`   (regenerated each run)
- `reports/design/kinetics-v0.5-review.pdf`   (Quarto-rendered)

The QMD is committed-friendly (deterministic from the TOML) so reviewers
can diff against earlier renders; the PDF is the deliverable.
"""

from __future__ import annotations

import re
import subprocess
import sys
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple


REPO_ROOT = Path(__file__).resolve().parents[2]
TOML_PATH = REPO_ROOT / "reports" / "params" / "calcium-v0.5.toml"
QMD_PATH = REPO_ROOT / "reports" / "design" / "kinetics-v0.5-review.qmd"
PDF_PATH = QMD_PATH.with_suffix(".pdf")
BIB_PATH = REPO_ROOT / "reports" / "params" / "calcium-v0.5-references.bib"


# ── Top-level namespace → chapter title mapping ──────────────────────
# Determines which TOML top-level keys become H1 chapter headings.
# Unknown keys fall through to "Other".
CHAPTER_TITLES = {
	"resting":   "Resting concentrations and calibration targets",
	"pm":        "Plasma-membrane Ca²⁺ leak",
	"ip3r":      "IP3R (deYoung-Keizer / Li-Rinzel)",
	"serca":     "SERCA cycle",
	"pmca":      "PMCA (basal + CaM-activated)",
	"cam":       "Calmodulin Ca²⁺ binding",
	"p2x1":      "P2X1 ATP-gated cation channel",
	"soce":      "Store-operated Ca²⁺ entry (STIM1 / MWC / puncta / Orai)",
	"ncx":       "Na⁺/Ca²⁺ exchanger",
	"agonists":  "Agonist forcing curves",
	"gpcr":      "GPCR cascade (P2Y1 / PAR1 / PAR4 / Gαq)",
	"pi_cycle":  "PI cycle (PLCβ / metabolism)",
	"mito":      "Mitochondrial Ca²⁺ (MCU + NCLX)",
	"buffers":   "Ca²⁺ buffers (cytosolic + DTS luminal)",
	"references": "References",
}


# ── Section block extracted from the TOML walk ───────────────────────


@dataclass
class Block:
	"""One TOML section as it appears in the file."""
	comments: List[str] = field(default_factory=list)     # raw stripped lines, paragraph-separated by ""
	section: Optional[str] = None                          # e.g. "ip3r.k_dyk"
	kv: List[Tuple[str, str, str]] = field(default_factory=list)  # (key, value, inline_comment)


# ── TOML parsing (raw, comment-preserving) ───────────────────────────


_SECTION_RE = re.compile(r"^\[(.+?)\]\s*$")
_KV_RE = re.compile(r"^\s*(\w+)\s*=\s*([^#]+?)(?:\s+#\s*(.*?))?\s*$")


_DIVIDER_RE = re.compile(r"^─+\s*(.+?)\s*─+$")


def parse_toml_blocks(path: Path) -> List[Block]:
	"""Walk the file linearly, grouping each comment block with the
	*following* section header. Comments appearing BETWEEN keys of a
	single section (rare) are discarded — the rendered table can't
	represent them cleanly, and the per-key inline ``# ...`` annotation
	is the more useful per-row text.
	"""
	blocks: List[Block] = []
	pending_comments: List[str] = []
	current_section: Optional[str] = None
	current_pre: List[str] = []
	current_kv: List[Tuple[str, str, str]] = []

	for raw in path.read_text(encoding='utf-8').splitlines():
		stripped = raw.strip()
		m_sec = _SECTION_RE.match(stripped)

		if m_sec and not stripped.startswith("[["):
			if current_section is not None:
				blocks.append(Block(comments=current_pre, section=current_section, kv=current_kv))
			current_section = m_sec.group(1)
			current_pre = pending_comments
			pending_comments = []
			current_kv = []

		elif stripped.startswith("#"):
			pending_comments.append(stripped.lstrip("#").lstrip())

		elif not stripped:
			if pending_comments and pending_comments[-1] != "":
				pending_comments.append("")

		else:
			m_kv = _KV_RE.match(raw)
			if m_kv:
				current_kv.append(
					(m_kv.group(1), m_kv.group(2).strip(), (m_kv.group(3) or "").strip())
				)
				pending_comments = []

	if current_section is not None:
		blocks.append(Block(comments=current_pre, section=current_section, kv=current_kv))

	return blocks


def extract_section_title(comments: List[str]) -> Optional[str]:
	"""If a leading comment line is a ``── X ──`` divider, return X."""
	for line in comments:
		if not line:
			continue
		m = _DIVIDER_RE.match(line)
		if m:
			return m.group(1)
		return None
	return None


# ── Reference resolution ─────────────────────────────────────────────


def load_references(path: Path) -> dict:
	"""Read [references.*] sub-tables and return a flat lookup dict."""
	with open(path, "rb") as f:
		data = tomllib.load(f)
	return data.get("references", {})


def build_match_index(refs: dict) -> List[Tuple[re.Pattern, str]]:
	"""Compile a regex index of citation strings → reference key.

	Each `[references.<key>]` may define a `match` list of strings that
	should hyperlink to it. We compile each into a word-boundary regex
	and return (pattern, ref_key) pairs in longest-match-first order
	so e.g. "Vu et al. 1991" wins over "Vu 1991".
	"""
	pairs: List[Tuple[str, str]] = []
	for key, meta in refs.items():
		for s in meta.get("match", []):
			pairs.append((s, key))
	# Longest first
	pairs.sort(key=lambda p: -len(p[0]))
	# Negative lookbehind also excludes `[` so a match inside an
	# already-substituted Markdown link `[Vu et al. 1991](#ref-...)` is
	# not re-matched on a subsequent pass — otherwise a shorter alias
	# like "Vu" or "Mahaut-Smith" would nest brackets inside the
	# longer-form link.
	return [(re.compile(r"(?<![\w\[])" + re.escape(s) + r"(?![\w])"), k) for s, k in pairs]


def ref_url(meta: dict) -> str:
	"""Resolve a reference's primary clickable URL."""
	doi = meta.get("doi")
	if doi:
		return f"https://doi.org/{doi}"
	return meta.get("url", "")


# ── Linkification ────────────────────────────────────────────────────


def linkify(text: str, match_index: List[Tuple[re.Pattern, str]], refs: dict) -> str:
	"""Convert citations / DOIs / URLs into Markdown hyperlinks."""

	# 1. Author-YEAR citations → anchor link into the [References]
	# section. The link target is always the in-document anchor; the
	# DOI/URL belongs in the bibliography entry itself, where readers
	# arrive after clicking the citation. We link whenever a matching
	# `[references.<key>]` exists, regardless of whether that entry has
	# a `doi` or `url` field set.
	def cite_sub(m, ref_key):
		txt = m.group(0)
		if ref_key in refs:
			return f"[{txt}](#ref-{ref_key})"
		return txt

	for pattern, key in match_index:
		text = pattern.sub(lambda m, k=key: cite_sub(m, k), text)

	# 2. Bare `doi:X` markers → clickable hyperlinks.
	text = re.sub(
		r"doi:\s*(10\.[^\s)\]\"']+)",
		r"[doi:\1](https://doi.org/\1)",
		text,
	)

	# 3. Bare URLs.
	text = re.sub(
		r"(?<![(\[\"'])(https?://[^\s)\]\"']+)",
		r"<\1>",  # angle-bracketed = clickable in pandoc/Quarto
		text,
	)

	return text


# ── Rendering ────────────────────────────────────────────────────────


def chapter_for(section: str) -> Tuple[str, str]:
	"""Return (top_level_key, chapter_title) for a section path."""
	top = section.split(".", 1)[0]
	return top, CHAPTER_TITLES.get(top, top.replace("_", " ").title())


def render_prose(comments: List[str], linkifier) -> str:
	"""Render a comment block as a sequence of paragraphs.

	``── X ──`` divider lines are dropped (their content is hoisted into
	the surrounding section's H2 by ``build_qmd``).
	"""
	# Filter divider lines out entirely.
	filtered = [line for line in comments if not _DIVIDER_RE.match(line)]
	# Collapse leading/trailing blanks.
	while filtered and filtered[-1] == "":
		filtered = filtered[:-1]
	while filtered and filtered[0] == "":
		filtered = filtered[1:]
	paragraphs: List[List[str]] = [[]]
	for line in filtered:
		if line == "":
			if paragraphs[-1]:
				paragraphs.append([])
		else:
			paragraphs[-1].append(line)
	out_paras = []
	for para in paragraphs:
		if not para:
			continue
		joined = " ".join(para)
		out_paras.append(linkifier(joined))
	return "\n\n".join(out_paras)


def render_kv_table(kv: List[Tuple[str, str, str]], linkifier) -> str:
	if not kv:
		return ""
	out = ["", "| Parameter | Value | Notes |", "|---|---|---|"]
	for key, val, comment in kv:
		notes = linkifier(comment) if comment else ""
		out.append(f"| `{key}` | `{val}` | {notes} |")
	out.append("")
	return "\n".join(out)


def _via_annotation(via_key: str, refs: dict) -> str:
	"""Format a `via = "<key>"` field as a clickable citation chain link.

	Returns ``" — *via* [Author YEAR](#ref-key)"`` when the upstream
	reference exists; falls back to a plain ``" — *via* <key>"``
	otherwise. The intent is for a reviewer scanning the bibliography
	to see at a glance "I learned about Bezprozvanny 1991 because
	Hoover 2011 cited it" and to click through to the upstream entry.
	"""
	upstream = refs.get(via_key)
	if upstream is None:
		return f" — *via* `{via_key}` (key not found in [references.*])"
	up_authors = upstream.get("authors", "")
	up_year = upstream.get("year", "")
	# Use a short "FirstAuthorSurname YEAR" form for the visible label.
	first_surname = up_authors.split(",")[0].split(" and ")[0].strip()
	label = f"{first_surname} {up_year}".strip()
	return f" — *via* [{label}](#ref-{via_key})"


def render_references(refs: dict) -> str:
	"""Render the bibliography section."""
	out = ["", "# References", ""]
	# Sort by first-author surname (best-effort; falls back to key).
	def sort_key(item):
		_, meta = item
		authors = meta.get("authors", "")
		return (authors.split(",")[0].strip(), meta.get("year", 0))
	for key, meta in sorted(refs.items(), key=sort_key):
		anchor = f"{{#ref-{key}}}"
		authors = meta.get("authors", "")
		year = meta.get("year", "n.d.")
		title = meta.get("title", "")
		journal = meta.get("journal", "")
		doi = meta.get("doi")
		url = meta.get("url")
		via = meta.get("via")
		link = (
			f"[doi:{doi}](https://doi.org/{doi})" if doi
			else (f"<{url}>" if url else "")
		)
		via_part = _via_annotation(via, refs) if via else ""
		entry = f"**{authors} ({year})** {anchor}. *{title}*. {journal}. {link}{via_part}"
		out.append("- " + entry)
	out.append("")
	return "\n".join(out)


# ── Top-level builder ────────────────────────────────────────────────


YAML_HEADER = """---
title: "Platelet WCM — calcium kinetics review (v0.5)"
subtitle: "Source-of-truth review document, auto-generated from `reports/params/calcium-v0.5.toml`"
author: "Steve Haigh"
date: today
format:
  pdf:
    pdf-engine: xelatex
    geometry: "margin=2.5cm"
    fontsize: 10pt
    mainfont: "STIX Two Text"
    monofont: "Menlo"
    colorlinks: true
    linkcolor: blue
    urlcolor: blue
    toc: true
    toc-depth: 3
    include-in-header: ../pandoc-header.tex
---

::: callout-note
This document is auto-generated from
[`reports/params/calcium-v0.5.toml`](../params/calcium-v0.5.toml).
Edit the TOML, then re-run
`runscripts/manual/buildKineticsReview.py` to regenerate.

Citations of the form *Author YEAR* are auto-linked to the
[References](#references) section if a matching entry exists in the
TOML's `[references.*]` block. Citations without a reference entry
appear as plain text — add a `[references.<key>]` block to the TOML
to make them clickable.
:::

"""


def build_qmd(blocks: List[Block], refs: dict) -> str:
	match_index = build_match_index(refs)
	def linker(text):
		return linkify(text, match_index, refs)

	out = [YAML_HEADER]

	current_chapter: Optional[str] = None
	for block in blocks:
		if block.section is None:
			continue
		if block.section == "references" or block.section.startswith("references."):
			continue  # Bibliography handled separately at the end.
		top, chapter_title = chapter_for(block.section)
		if top != current_chapter:
			current_chapter = top
			out.append(f"\n# {chapter_title}\n")
		# Use the ── X ── divider title (if present in the leading
		# comments) as the H2 subtitle for the section.
		section_title = extract_section_title(block.comments)
		if section_title:
			out.append(f"## `[{block.section}]` — {section_title}\n")
		else:
			out.append(f"## `[{block.section}]`\n")
		if block.comments:
			out.append(render_prose(block.comments, linker))
		out.append(render_kv_table(block.kv, linker))

	if refs:
		out.append(render_references(refs))

	return "\n".join(out)


def write_bibtex(refs: dict, path: Path) -> None:
	"""Emit a BibTeX file from the parsed `[references.*]` entries.

	The output is import-friendly for Zotero (drag-and-drop) and is
	also re-usable as a `\\bibliography{...}` source for any future
	LaTeX writeup. One `@article{key, ...}` per reference key.
	"""
	lines = [
		"% Auto-generated from reports/params/calcium-v0.5.toml [references.*].",
		"% Regenerate with: PYTHONPATH=$PWD python runscripts/manual/buildKineticsReview.py",
		"",
	]
	def escape_bibtex(s: str) -> str:
		"""Escape characters that would break a BibTeX `{...}` value.

		Bare `{` / `}` in a value collide with BibTeX's brace grouping;
		`\\` would be interpreted as a TeX control sequence. Prefix each
		with a backslash. No field in the current TOML triggers any of
		these, but the escaping is cheap insurance against future
		additions (e.g., titles or notes containing math braces).
		"""
		return s.replace("\\", r"\\").replace("{", r"\{").replace("}", r"\}")

	for key in sorted(refs.keys()):
		meta = refs[key]
		fields: List[Tuple[str, str]] = []
		for f in ("authors", "title", "journal", "year", "volume", "pages", "doi", "url"):
			v = meta.get(f)
			if v is None or v == "":
				continue
			bibtex_field = {"authors": "author"}.get(f, f)
			fields.append((bibtex_field, escape_bibtex(str(v))))
		# Record citation route in a BibTeX `note` field so the chain
		# survives Zotero import / future LaTeX writeups.
		via = meta.get("via")
		if via:
			upstream = refs.get(via, {})
			up_authors = upstream.get("authors", "")
			up_year = upstream.get("year", "")
			first_surname = up_authors.split(",")[0].split(" and ")[0].strip()
			fields.append(("note", escape_bibtex(f"via {first_surname} {up_year}".strip())))
		lines.append(f"@article{{{key},")
		for k, v in fields:
			lines.append(f"  {k:8s} = {{{v}}},")
		lines.append("}")
		lines.append("")
	path.write_text("\n".join(lines), encoding='utf-8')


def main() -> int:
	if not TOML_PATH.exists():
		print(f"TOML not found at {TOML_PATH}", file=sys.stderr)
		return 1

	refs = load_references(TOML_PATH)
	blocks = parse_toml_blocks(TOML_PATH)
	qmd_text = build_qmd(blocks, refs)

	QMD_PATH.parent.mkdir(parents=True, exist_ok=True)
	QMD_PATH.write_text(qmd_text, encoding='utf-8')
	print(f"Wrote {QMD_PATH}")

	# BibTeX side-output for Zotero import / future LaTeX use.
	if refs:
		write_bibtex(refs, BIB_PATH)
		print(f"Wrote {BIB_PATH} ({len(refs)} reference{'s' if len(refs) != 1 else ''})")

	# Render via Quarto. The render command will overwrite the PDF in
	# place. We capture stderr only on failure so a successful run is
	# quiet.
	result = subprocess.run(
		["quarto", "render", str(QMD_PATH), "--to", "pdf"],
		cwd=REPO_ROOT,
		capture_output=True,
		text=True,
	)
	if result.returncode != 0:
		print("Quarto render failed:", file=sys.stderr)
		print(result.stdout, file=sys.stderr)
		print(result.stderr, file=sys.stderr)
		return result.returncode

	print(f"Rendered {PDF_PATH}")
	return 0


if __name__ == "__main__":
	sys.exit(main())
