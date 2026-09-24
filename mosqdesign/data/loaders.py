"""FASTA loaders (files fetched live by scripts/fetch_data.sh)."""
from __future__ import annotations

import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data")


def read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs
