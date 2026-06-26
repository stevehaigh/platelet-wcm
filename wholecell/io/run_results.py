"""Locate variant + cell output directories within a simulation run.

These helpers walk the wcEcoli-inherited output nesting
(``{variant}/{seed}/generation_{n}/{daughter}/simOut``) so analysis code can
discover what a run produced. They were previously part of the (now-removed)
Dash webapp's ``results`` module; ``analysisPlatelet.py`` is the remaining
consumer.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List

VARIANT_PATTERN = re.compile(r'.+_\d{6}')
SEED_PATTERN = re.compile(r'\d{6}')
GENERATION_PATTERN = re.compile(r'generation_\d{6}')
DAUGHTER_PATTERN = re.compile(r'\d{6}')


def find_variants(sim_dir: str) -> List[str]:
	"""Find variant subdirectories within a simulation output directory."""

	variants = []
	if not os.path.isdir(sim_dir):
		return variants

	for name in sorted(os.listdir(sim_dir)):
		path = os.path.join(sim_dir, name)
		if os.path.isdir(path) and VARIANT_PATTERN.fullmatch(name):
			variants.append(name)
	return variants


def find_cells(sim_dir: str, variant: str) -> List[Dict[str, str]]:
	"""Find all cell simulation directories for a given variant.

	Returns a list of dicts with keys: seed, generation, daughter, simout_path.
	"""

	cells = []
	variant_dir = os.path.join(sim_dir, variant)
	if not os.path.isdir(variant_dir):
		return cells

	for seed in sorted(os.listdir(variant_dir)):
		if not SEED_PATTERN.fullmatch(seed):
			continue
		seed_dir = os.path.join(variant_dir, seed)
		if not os.path.isdir(seed_dir):
			continue

		for gen in sorted(os.listdir(seed_dir)):
			if not GENERATION_PATTERN.fullmatch(gen):
				continue
			gen_dir = os.path.join(seed_dir, gen)
			if not os.path.isdir(gen_dir):
				continue

			for daughter in sorted(os.listdir(gen_dir)):
				if not DAUGHTER_PATTERN.fullmatch(daughter):
					continue
				simout = os.path.join(gen_dir, daughter, 'simOut')
				if os.path.isdir(simout):
					cells.append({
						'seed': seed,
						'generation': gen,
						'daughter': daughter,
						'simout_path': simout,
					})
	return cells
