"""File helpers shared by recognition, evidence, and export."""

import csv
import json
from hashlib import file_digest
from pathlib import Path


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def save(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":")), encoding="utf-8")
    temp.replace(path)


def sha(path):
    with Path(path).open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest().upper()


def write_csv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def checked(path, digest):
    path = Path(path).resolve(strict=True)
    if sha(path) != digest:
        raise ValueError(f"Fingerprint mismatch: {path}")
    return path
