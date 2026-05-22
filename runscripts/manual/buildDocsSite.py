"""Render every doc under reports/ to HTML via Quarto + build an index page.

Quarto handles both .md and .qmd files. Each source file under
reports/{data,decks,design,external,lab-books} (plus any reports/*.md at
the top level) is rendered to HTML into reports/site/<category>/, and a
reports/site/index.html is emitted listing every rendered doc grouped by
category with title, source path, and a link.

The site directory is gitignored by default; flip the `/reports/site/`
line in .gitignore (and `git add reports/site/`) to publish via the
already-configured GitHub Pages instance at
`https://stevehaigh.github.io/platelet-wcm/reports/site/index.html`.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/buildDocsSite.py
    PYTHONPATH=$PWD python runscripts/manual/buildDocsSite.py --only design
    PYTHONPATH=$PWD python runscripts/manual/buildDocsSite.py --jobs 4
"""

from __future__ import annotations

import argparse
import concurrent.futures
import html
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REPORTS = REPO_ROOT / 'reports'
SITE = REPORTS / 'site'

# Order controls how categories appear in the index.
CATEGORIES = ['design', 'lab-books', 'data', 'external', 'decks']
ROOT_LABEL = 'top-level'

DATE_RE = re.compile(r'(20\d{2}-\d{2}-\d{2})')
# Minimal frontmatter key-value lines: `key: value` or `key: "value"` at column 0.
_FRONT_KV_RE = re.compile(r'^([A-Za-z_]+)\s*:\s*(.+?)\s*$')


def _read_frontmatter(text: str) -> dict[str, str]:
	"""Extract top-level `key: value` lines from a YAML frontmatter block.

	Deliberately minimal — only handles the flat string-valued keys we
	care about (`title`, `date`, `author`). Avoids a PyYAML dependency.
	"""
	if not text.startswith('---'):
		return {}
	parts = text.split('---', 2)
	if len(parts) < 3:
		return {}
	out: dict[str, str] = {}
	for line in parts[1].splitlines():
		if not line or line.lstrip().startswith('#'):
			continue
		m = _FRONT_KV_RE.match(line)
		if not m:
			continue
		key, value = m.group(1), m.group(2)
		# Strip surrounding quotes if any.
		if (value.startswith('"') and value.endswith('"')) or (
				value.startswith("'") and value.endswith("'")):
			value = value[1:-1]
		out[key] = value
	return out


@dataclass
class Doc:
	category: str       # one of CATEGORIES or ROOT_LABEL
	src: Path           # absolute path to .md / .qmd
	html_rel: Path      # path relative to SITE/, set after render
	title: str          # from frontmatter or first H1
	date: str           # YYYY-MM-DD if extractable, else ''


# ── Discovery ─────────────────────────────────────────────────────────────

def find_docs(only: list[str] | None = None) -> list[tuple[str, Path]]:
	"""Return (category, source_path) pairs for every doc to render."""
	pairs: list[tuple[str, Path]] = []
	cats_to_scan = only or [ROOT_LABEL] + CATEGORIES
	# Top-level reports/*.md and *.qmd (e.g. dissertation-notes.md).
	if ROOT_LABEL in cats_to_scan:
		for p in sorted(REPORTS.glob('*.md')) + sorted(REPORTS.glob('*.qmd')):
			if p.is_file():
				pairs.append((ROOT_LABEL, p))
	# Categorised subdirs.
	for cat in CATEGORIES:
		if cat not in cats_to_scan:
			continue
		cat_dir = REPORTS / cat
		if not cat_dir.is_dir():
			continue
		for p in sorted(cat_dir.glob('*.md')) + sorted(cat_dir.glob('*.qmd')):
			if p.is_file():
				pairs.append((cat, p))
	return pairs


# ── Metadata extraction ───────────────────────────────────────────────────

def extract_title(path: Path) -> str:
	"""Frontmatter title, falling back to first H1, then filename."""
	try:
		text = path.read_text(encoding='utf-8')
	except (OSError, UnicodeDecodeError):
		return path.stem
	fm = _read_frontmatter(text)
	if fm.get('title'):
		return fm['title']
	for line in text.splitlines():
		if line.startswith('# '):
			return line[2:].strip()
	return path.stem


def extract_date(path: Path) -> str:
	"""YYYY-MM-DD from filename or frontmatter; '' if neither has one."""
	m = DATE_RE.search(path.stem)
	if m:
		return m.group(1)
	try:
		text = path.read_text(encoding='utf-8')
	except (OSError, UnicodeDecodeError):
		return ''
	fm = _read_frontmatter(text)
	return fm.get('date', '')


# ── Render ────────────────────────────────────────────────────────────────

def render(src: Path, top_out_dir: Path, log_lines: list[str]) -> bool:
	"""quarto render <src> --to html --output-dir <top_out_dir>.

	The repo's root `_quarto.yml` makes Quarto treat the whole repo as
	a project, which means `--output-dir X` produces output at
	`X/<src.relative_to(repo_root)>.html` regardless of `cwd`. We
	accept that layout here and flatten it post-render via
	`_flatten_quarto_output()` so the user-visible path is the
	expected `<top_out_dir>/<category>/<stem>.html`.
	"""
	top_out_dir.mkdir(parents=True, exist_ok=True)
	cmd = ['quarto', 'render', str(src.relative_to(REPO_ROOT)),
		'--to', 'html', '--output-dir', str(top_out_dir.resolve())]
	start = time.time()
	result = subprocess.run(cmd, capture_output=True, text=True,
		cwd=str(REPO_ROOT))
	elapsed = time.time() - start
	rel = src.relative_to(REPO_ROOT)
	if result.returncode != 0:
		log_lines.append(
			f'  ✗ {rel} ({elapsed:.1f} s)\n    {result.stderr.strip()[-400:]}')
		return False
	log_lines.append(f'  ✓ {rel} ({elapsed:.1f} s)')
	return True


def _flatten_quarto_output(out_dir: Path) -> None:
	"""Lift `<out_dir>/reports/*` up to `<out_dir>/*`.

	Quarto preserves the source's project-relative path inside the
	output dir, so renders of `reports/lab-books/foo.md` land at
	`<out_dir>/reports/lab-books/foo.html`. We expose them at
	`<out_dir>/lab-books/foo.html` for sanity.
	"""
	inner = out_dir / 'reports'
	if not inner.is_dir():
		return
	for entry in inner.iterdir():
		dest = out_dir / entry.name
		if dest.exists() and dest.is_dir() and entry.is_dir():
			# Merge sibling dir into existing (e.g. category collisions).
			for child in entry.iterdir():
				shutil.move(str(child), str(dest / child.name))
			entry.rmdir()
		else:
			shutil.move(str(entry), str(dest))
	# `inner` should now be empty.
	try:
		inner.rmdir()
	except OSError:
		# Best-effort: if something else landed here, leave it.
		pass


# ── Index page ────────────────────────────────────────────────────────────

CATEGORY_LABELS = {
	ROOT_LABEL: 'Top-level',
	'design':    'Design docs',
	'lab-books': 'Lab books',
	'data':      'Data provenance',
	'external':  'External references',
	'decks':     'Slide decks',
}

INDEX_CSS = """
* { box-sizing: border-box; }
body {
	font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
	margin: 0; padding: 40px 20px; color: #1f2328;
	background: #ffffff; line-height: 1.5;
}
.container { max-width: 1080px; margin: 0 auto; }
h1 { font-size: 28px; margin: 0 0 4px 0; }
.subtitle { color: #57606a; font-size: 14px; margin-bottom: 32px; }
h2 {
	font-size: 18px; margin: 36px 0 8px 0;
	padding-bottom: 6px; border-bottom: 1px solid #d0d7de;
}
.cat-count { color: #57606a; font-weight: 400; font-size: 13px; }
table { border-collapse: collapse; width: 100%; font-size: 14px; }
th, td { text-align: left; padding: 8px 12px; vertical-align: top; }
th {
	color: #57606a; font-weight: 500; font-size: 12px;
	text-transform: uppercase; letter-spacing: 0.04em;
	border-bottom: 1px solid #d0d7de;
}
tr:nth-child(even) { background: #f6f8fa; }
td.date { color: #57606a; white-space: nowrap; width: 100px; font-variant-numeric: tabular-nums; }
td.src { color: #57606a; font-family: ui-monospace, "SF Mono", Menlo, monospace; font-size: 12px; }
a { color: #0969da; text-decoration: none; }
a:hover { text-decoration: underline; }
footer {
	margin-top: 48px; padding-top: 20px;
	border-top: 1px solid #d0d7de;
	color: #57606a; font-size: 13px;
}
"""


def build_index(docs: list[Doc], index_path: Path, generated_at: str) -> None:
	"""Emit reports/site/index.html with one table per category."""
	by_cat: dict[str, list[Doc]] = defaultdict(list)
	for d in docs:
		by_cat[d.category].append(d)

	# Sort: lab books newest-first by date, everything else by title.
	for cat, items in by_cat.items():
		if cat == 'lab-books':
			items.sort(key=lambda d: (d.date or '', d.title), reverse=True)
		else:
			items.sort(key=lambda d: (d.date or '', d.title), reverse=True)

	cat_order = [ROOT_LABEL] + CATEGORIES
	sections: list[str] = []
	for cat in cat_order:
		if cat not in by_cat:
			continue
		items = by_cat[cat]
		label = CATEGORY_LABELS.get(cat, cat)
		rows = []
		for d in items:
			date_cell = html.escape(d.date) if d.date else ''
			title_cell = (
				f'<a href="{html.escape(str(d.html_rel))}">'
				f'{html.escape(d.title)}</a>'
			)
			src_cell = html.escape(str(d.src.relative_to(REPORTS)))
			rows.append(
				f'<tr><td class="date">{date_cell}</td>'
				f'<td>{title_cell}</td>'
				f'<td class="src">{src_cell}</td></tr>')
		sections.append(
			f'<h2>{html.escape(label)} '
			f'<span class="cat-count">({len(items)})</span></h2>\n'
			f'<table>\n'
			f'<thead><tr><th>Date</th><th>Title</th><th>Source</th></tr></thead>\n'
			f'<tbody>\n{chr(10).join(rows)}\n</tbody>\n</table>')

	page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Platelet WCM — docs index</title>
<style>{INDEX_CSS}</style>
</head>
<body>
<div class="container">
<h1>Platelet WCM — docs index</h1>
<p class="subtitle">{len(docs)} documents · generated {html.escape(generated_at)} ·
<a href="https://github.com/stevehaigh/platelet-wcm">source</a></p>
{chr(10).join(sections)}
<footer>Built by <code>runscripts/manual/buildDocsSite.py</code>. Re-run
to refresh; the index, the per-doc HTML, and all supporting assets
live under <code>reports/site/</code>.</footer>
</div>
</body>
</html>
"""
	index_path.write_text(page, encoding='utf-8')


# ── CLI ───────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--out-dir', default=str(SITE),
		help=f'Output directory. Default: {SITE.relative_to(REPO_ROOT)}/')
	parser.add_argument('--only', action='append', default=None,
		help='Restrict to one category (design / lab-books / data / external / '
		'decks / top-level). May be repeated.')
	parser.add_argument('--jobs', type=int, default=1,
		help='Parallel quarto renders. Default 1 (Quarto can be flaky in '
		'parallel for some projects).')
	args = parser.parse_args(argv)

	out_dir = Path(args.out_dir)
	if not out_dir.is_absolute():
		out_dir = REPO_ROOT / out_dir

	pairs = find_docs(only=args.only)
	if not pairs:
		print('No docs found.', file=sys.stderr)
		return 1
	print(f'Rendering {len(pairs)} doc(s) into {out_dir.relative_to(REPO_ROOT)}/ …\n')

	results: list[Doc] = []
	log_lines: list[str] = []

	def _do(pair: tuple[str, Path]) -> tuple[bool, str, Path]:
		cat, src = pair
		# All renders share the same top-level output dir; Quarto's
		# project detection preserves the source's relative path
		# inside that dir, and we flatten it post-render.
		ok = render(src, out_dir, log_lines)
		return ok, cat, src

	if args.jobs > 1:
		with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
			outcomes = list(ex.map(_do, pairs))
	else:
		outcomes = [_do(p) for p in pairs]

	for line in log_lines:
		print(line)

	# Flatten Quarto's nested layout once, after every render is done.
	_flatten_quarto_output(out_dir)

	for ok, cat, src in outcomes:
		if not ok:
			continue
		# After flattening: <out_dir>/<src.parent.name>/<stem>.html for
		# categorised docs, <out_dir>/<stem>.html for top-level.
		if cat == ROOT_LABEL:
			html_path = out_dir / (src.stem + '.html')
		else:
			html_path = out_dir / cat / (src.stem + '.html')
		if not html_path.exists():
			# Defensive: glob for the file in case Quarto renamed it.
			candidates = list(html_path.parent.glob(f'{src.stem}*.html'))
			if candidates:
				html_path = candidates[0]
		results.append(Doc(
			category=cat, src=src,
			html_rel=html_path.relative_to(out_dir),
			title=extract_title(src),
			date=extract_date(src),
		))

	index_path = out_dir / 'index.html'
	build_index(results, index_path,
		generated_at=datetime.now().strftime('%Y-%m-%d %H:%M'))

	n_ok = len(results)
	n_failed = sum(1 for ok, _, _ in outcomes if not ok)
	print(f'\nRendered {n_ok}; failed {n_failed}.')
	print(f'Index: {index_path}')
	if n_failed:
		print(f'(See the ✗ lines above for failure details.)')
	return 0 if n_failed == 0 else 2


if __name__ == '__main__':
	sys.exit(main())
