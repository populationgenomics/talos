"""
SQLite database layer for Talos web application.
Stores analyst decisions and run registrations.
"""

import sqlite3
from contextlib import contextmanager
from pathlib import Path

from loguru import logger

SCHEMA_VERSION = 1

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS runs (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    label       TEXT NOT NULL,
    dataset     TEXT NOT NULL DEFAULT '',
    json_path   TEXT NOT NULL UNIQUE,
    run_date    TEXT NOT NULL DEFAULT '',
    registered  TEXT NOT NULL DEFAULT (datetime('now')),
    notes       TEXT NOT NULL DEFAULT ''
);

CREATE TABLE IF NOT EXISTS decisions (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sample          TEXT NOT NULL,
    chrom           TEXT NOT NULL,
    pos             INTEGER NOT NULL,
    ref             TEXT NOT NULL,
    alt             TEXT NOT NULL,
    gene            TEXT NOT NULL,
    status          TEXT NOT NULL DEFAULT 'pending',
    acmg_class      TEXT,
    analyst         TEXT NOT NULL DEFAULT '',
    comment         TEXT NOT NULL DEFAULT '',
    created_at      TEXT NOT NULL DEFAULT (datetime('now')),
    updated_at      TEXT NOT NULL DEFAULT (datetime('now')),
    UNIQUE(sample, chrom, pos, ref, alt, gene)
);

CREATE TABLE IF NOT EXISTS decision_history (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    decision_id     INTEGER NOT NULL REFERENCES decisions(id),
    status          TEXT NOT NULL,
    acmg_class      TEXT,
    analyst         TEXT NOT NULL DEFAULT '',
    comment         TEXT NOT NULL DEFAULT '',
    created_at      TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX IF NOT EXISTS idx_decisions_variant ON decisions(sample, chrom, pos, ref, alt, gene);
CREATE INDEX IF NOT EXISTS idx_decisions_status ON decisions(status);
CREATE INDEX IF NOT EXISTS idx_history_decision ON decision_history(decision_id);
"""


def init_db(db_path: str) -> sqlite3.Connection:
    """Initialise the database, creating tables if needed. Returns a connection in WAL mode."""
    path = Path(db_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(str(path), check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute('PRAGMA journal_mode=WAL')
    conn.execute('PRAGMA foreign_keys=ON')
    conn.executescript(SCHEMA_SQL)
    conn.commit()
    logger.info(f'Database initialised at {path}')
    return conn


@contextmanager
def transaction(conn: sqlite3.Connection):
    """Context manager for a database transaction."""
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise


# --- Run CRUD ---


def list_runs(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    return conn.execute('SELECT * FROM runs ORDER BY registered DESC').fetchall()


def get_run(conn: sqlite3.Connection, run_id: int) -> sqlite3.Row | None:
    return conn.execute('SELECT * FROM runs WHERE id = ?', (run_id,)).fetchone()


def register_run(
    conn: sqlite3.Connection,
    label: str,
    json_path: str,
    dataset: str = '',
    run_date: str = '',
    notes: str = '',
) -> int:
    with transaction(conn):
        cursor = conn.execute(
            'INSERT INTO runs (label, json_path, dataset, run_date, notes) VALUES (?, ?, ?, ?, ?)',
            (label, json_path, dataset, run_date, notes),
        )
        return cursor.lastrowid


def delete_run(conn: sqlite3.Connection, run_id: int) -> bool:
    with transaction(conn):
        cursor = conn.execute('DELETE FROM runs WHERE id = ?', (run_id,))
        return cursor.rowcount > 0


# --- Decision CRUD ---


VALID_STATUSES = frozenset({'pending', 'under_review', 'reportable', 'uncertain', 'not_reportable', 'solved'})
VALID_ACMG = frozenset({None, 'pathogenic', 'likely_pathogenic', 'uncertain_significance', 'likely_benign', 'benign'})


def upsert_decision(
    conn: sqlite3.Connection,
    *,
    sample: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene: str,
    status: str,
    analyst: str = '',
    comment: str = '',
    acmg_class: str | None = None,
) -> dict:
    """Create or update a decision. Appends to decision_history on every change."""
    if status not in VALID_STATUSES:
        raise ValueError(f'Invalid status: {status}. Must be one of {VALID_STATUSES}')
    if acmg_class not in VALID_ACMG:
        raise ValueError(f'Invalid ACMG class: {acmg_class}. Must be one of {VALID_ACMG}')

    with transaction(conn):
        # upsert the decision
        conn.execute(
            """
            INSERT INTO decisions (sample, chrom, pos, ref, alt, gene, status, acmg_class, analyst, comment)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(sample, chrom, pos, ref, alt, gene) DO UPDATE SET
                status = excluded.status,
                acmg_class = excluded.acmg_class,
                analyst = excluded.analyst,
                comment = excluded.comment,
                updated_at = datetime('now')
            """,
            (sample, chrom, pos, ref, alt, gene, status, acmg_class, analyst, comment),
        )

        # get the decision id
        row = conn.execute(
            'SELECT id FROM decisions WHERE sample=? AND chrom=? AND pos=? AND ref=? AND alt=? AND gene=?',
            (sample, chrom, pos, ref, alt, gene),
        ).fetchone()

        # append to history
        conn.execute(
            'INSERT INTO decision_history (decision_id, status, acmg_class, analyst, comment) VALUES (?, ?, ?, ?, ?)',
            (row['id'], status, acmg_class, analyst, comment),
        )

    return dict(row)


def get_decision(
    conn: sqlite3.Connection,
    *,
    sample: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene: str,
) -> dict | None:
    row = conn.execute(
        'SELECT * FROM decisions WHERE sample=? AND chrom=? AND pos=? AND ref=? AND alt=? AND gene=?',
        (sample, chrom, pos, ref, alt, gene),
    ).fetchone()
    return dict(row) if row else None


def get_decisions_for_keys(conn: sqlite3.Connection, keys: list[tuple]) -> dict[tuple, dict]:
    """Batch lookup decisions by natural keys. Returns a dict of key -> decision row."""
    if not keys:
        return {}

    result = {}
    # SQLite has a variable limit, so batch in chunks
    batch_size = 100
    for i in range(0, len(keys), batch_size):
        batch = keys[i : i + batch_size]
        placeholders = ' OR '.join(['(sample=? AND chrom=? AND pos=? AND ref=? AND alt=? AND gene=?)'] * len(batch))
        params = [p for key in batch for p in key]
        rows = conn.execute(f'SELECT * FROM decisions WHERE {placeholders}', params).fetchall()  # noqa: S608
        for row in rows:
            key = (row['sample'], row['chrom'], row['pos'], row['ref'], row['alt'], row['gene'])
            result[key] = dict(row)
    return result


def get_decision_history(conn: sqlite3.Connection, decision_id: int) -> list[dict]:
    rows = conn.execute(
        'SELECT * FROM decision_history WHERE decision_id = ? ORDER BY created_at ASC',
        (decision_id,),
    ).fetchall()
    return [dict(r) for r in rows]


def get_decision_counts(conn: sqlite3.Connection) -> dict[str, int]:
    """Get total counts of decisions by status."""
    rows = conn.execute('SELECT status, COUNT(*) as count FROM decisions GROUP BY status').fetchall()
    return {row['status']: row['count'] for row in rows}
