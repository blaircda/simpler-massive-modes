"""
Shared setup for the uplift-structure explorations.

These scripts only *read* the computational output of so10_uplift/
(data/R10.json, data/r9.json) and reuse its library modules (uplift.py,
weights.py).  They never import baseline files and never modify the data.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SO10 = os.path.join(ROOT, "so10_uplift")
DATA = os.path.join(SO10, "data")
RESULTS = os.path.join(HERE, "results")
if SO10 not in sys.path:
    sys.path.insert(0, SO10)


def load_R10():
    with open(os.path.join(DATA, "R10.json")) as f:
        return json.load(f)["levels"]


def as_dict(entries):
    """JSON irrep list -> {doubled ambient weight tuple: multiplicity}."""
    return {tuple(e["ambient2"]): e["mult"] for e in entries}


def R9(levels, n):
    return as_dict(levels[str(n)]["R9_B4"])


def canonical(levels, n):
    return as_dict(levels[str(n)]["canonical_R10_D5"])


def extremal(levels, n):
    return {k: as_dict(v) for k, v in levels[str(n)]["extremal_R10_D5"].items()}


def max_level(levels):
    return max(int(k) for k in levels)
