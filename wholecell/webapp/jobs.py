"""Job manager — SQLite job queue and subprocess runner for platelet sims."""

from __future__ import annotations

import json
import os
import re
import sqlite3
import subprocess
import threading
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional


DB_FILENAME = 'webapp_jobs.db'
PHASES = ['queued', 'simulating', 'analyzing', 'done', 'failed']

# Typical durations for progress estimation (seconds)
PHASE_DURATIONS = {
	'simulating': 30,
	'analyzing': 10,
}


ALLOWED_COLUMNS = frozenset([
	'status', 'output_dir', 'pid', 'started_at', 'finished_at', 'error_message',
])


class JobManager:
	"""Simple job queue backed by SQLite."""

	def __init__(self, db_path: str, repo_root: str,
			docker_image: str = 'steve-wcm-code') -> None:
		self.db_path = db_path
		self.repo_root = repo_root
		self.docker_image = docker_image
		self._init_db()
		self._worker_thread: Optional[threading.Thread] = None
		self._worker_lock = threading.Lock()
		self._stop_event = threading.Event()

	def _init_db(self) -> None:
		with self._connect() as conn:
			conn.execute('''CREATE TABLE IF NOT EXISTS jobs (
				id INTEGER PRIMARY KEY AUTOINCREMENT,
				status TEXT DEFAULT 'queued',
				variant TEXT,
				config_json TEXT,
				description TEXT,
				pid INTEGER,
				output_dir TEXT,
				started_at TEXT,
				finished_at TEXT,
				error_message TEXT
			)''')

	def _connect(self) -> sqlite3.Connection:
		conn = sqlite3.connect(self.db_path, timeout=30)
		conn.row_factory = sqlite3.Row
		conn.execute('PRAGMA journal_mode=WAL')
		return conn

	def submit(self, config: Dict[str, Any]) -> int:
		"""Submit a new job. Returns the job ID."""

		now = datetime.now(timezone.utc).isoformat()
		with self._connect() as conn:
			cursor = conn.execute(
				'INSERT INTO jobs (status, variant, config_json, description, started_at) '
				'VALUES (?, ?, ?, ?, ?)',
				('queued', 'platelet',
				 json.dumps(config), config.get('description', ''), now))
			job_id = cursor.lastrowid

		self._ensure_worker()
		return job_id

	def list_jobs(self) -> List[Dict[str, Any]]:
		"""List all jobs, newest first."""

		with self._connect() as conn:
			rows = conn.execute(
				'SELECT * FROM jobs ORDER BY id DESC').fetchall()
		return [dict(row) for row in rows]

	def get_job(self, job_id: int) -> Optional[Dict[str, Any]]:
		with self._connect() as conn:
			row = conn.execute(
				'SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone()
		return dict(row) if row else None

	def _update_status(self, job_id: int, status: str, **kwargs) -> None:
		fields = ['status = ?']
		values = [status]
		for k, v in kwargs.items():
			if k not in ALLOWED_COLUMNS:
				raise ValueError(f'Invalid column name: {k}')
			fields.append(f'{k} = ?')
			values.append(v)
		values.append(job_id)

		with self._connect() as conn:
			conn.execute(
				f'UPDATE jobs SET {", ".join(fields)} WHERE id = ?', values)

	def _ensure_worker(self) -> None:
		with self._worker_lock:
			if self._worker_thread is None or not self._worker_thread.is_alive():
				self._stop_event.clear()
				self._worker_thread = threading.Thread(
					target=self._worker_loop, daemon=True)
				self._worker_thread.start()

	def _worker_loop(self) -> None:
		"""Process queued jobs one at a time, polling when idle."""

		while not self._stop_event.is_set():
			with self._connect() as conn:
				# Atomically claim the next queued job to prevent duplicates
				conn.execute(
					'UPDATE jobs SET status = ? WHERE id = '
					'(SELECT id FROM jobs WHERE status = ? ORDER BY id LIMIT 1)',
					('simulating', 'queued'))
				if conn.total_changes == 0:
					self._stop_event.wait(timeout=5)
					continue
				row = conn.execute(
					'SELECT * FROM jobs WHERE status = ? ORDER BY id LIMIT 1',
					('simulating',)).fetchone()

			job = dict(row)
			self._run_job(job)

	def _run_job(self, job: Dict[str, Any]) -> None:
		"""Execute a platelet simulation job — sim then analysis plots.

		When running inside a container (PLATELET_WEBAPP_MODE=container),
		uses local Python subprocess. Otherwise, shells out to Docker.
		"""

		job_id = job['id']
		config = json.loads(job['config_json'])
		timestamp = datetime.now().strftime('%Y%m%d.%H%M%S')
		desc = re.sub(r'[^a-zA-Z0-9_-]', '', config.get('description', '')) or 'webapp'
		sim_outdir = f'{timestamp}___{desc}'
		out_path = os.path.join(self.repo_root, 'out', sim_outdir)

		# Verify the resolved path stays inside the out/ directory
		out_root = os.path.realpath(os.path.join(self.repo_root, 'out'))
		if not os.path.realpath(out_path).startswith(out_root + os.sep):
			raise ValueError(f'Invalid output path: {out_path}')

		self._update_status(job_id, 'simulating', output_dir=out_path,
			pid=os.getpid())

		in_container = os.environ.get('PLATELET_WEBAPP_MODE') == 'container'
		length_sec = int(config.get('length_sec', 60))
		seed = int(config.get('seed', 0))

		try:
			# Phase 1: Simulation
			cmd_sim = self._build_cmd(
				[
					'python', 'runscripts/manual/runPlateletSim.py',
					sim_outdir,
					'--length', str(length_sec),
					'--seed', str(seed),
					'--no-log-to-shell',
				],
				in_container, out_root)
			proc = subprocess.run(cmd_sim, capture_output=True, text=True,
				cwd=self.repo_root if in_container else None)
			if proc.returncode != 0:
				raise RuntimeError(f'Platelet sim failed:\n{proc.stderr[-2000:]}')

			# Phase 2: Analysis
			self._update_status(job_id, 'analyzing')
			cmd_analysis = self._build_cmd(
				['python', 'runscripts/manual/analysisPlatelet.py', sim_outdir],
				in_container, out_root)
			proc = subprocess.run(cmd_analysis, capture_output=True, text=True,
				cwd=self.repo_root if in_container else None)
			if proc.returncode != 0:
				raise RuntimeError(f'Platelet analysis failed:\n{proc.stderr[-2000:]}')

			now = datetime.now(timezone.utc).isoformat()
			self._update_status(job_id, 'done', finished_at=now)
		except Exception as e:
			now = datetime.now(timezone.utc).isoformat()
			self._update_status(
				job_id, 'failed', finished_at=now,
				error_message=str(e)[:4000])

	def _build_cmd(self, args: List[str], in_container: bool,
			out_root: str) -> List[str]:
		"""Build a command, wrapping with docker run if not in a container."""
		if in_container:
			return args
		return [
			'docker', 'run', '--rm',
			'-v', f'{out_root}:/platelet-wcm/out',
			self.docker_image,
		] + args

	def stop(self) -> None:
		self._stop_event.set()
