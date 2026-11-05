"""Reading tab-separated tables and the numbers in them."""
import csv
import os


def read_tsv(path, missing_ok=False):
    """Rows of a tab-separated file with a header, as dicts; lines starting with "#" are skipped.
    missing_ok: a missing or empty file gives [] instead of an error."""
    if missing_ok and (not os.path.exists(path) or os.path.getsize(path) == 0):
        return []
    with open(path) as f:
        return list(csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t"))


def num(x):
    """A number, or None for empty, NA, NaN or text."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v
