#!/usr/bin/env python3
"""Fetch the small Eurostat LE65 validation slice as a reproducible CSV."""

import csv
import itertools
import math
import sys
from pathlib import Path

import requests


API = "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/demo_mlexpec"
GEOS = ["EU27_2020", "ES", "SE", "PL"]
YEARS = ["2019", "2023"]


def ordered_codes(dimension):
    index = dimension["category"]["index"]
    if isinstance(index, list):
        return index
    return [code for code, _ in sorted(index.items(), key=lambda item: item[1])]


def flatten_jsonstat(payload):
    ids = payload["id"]
    sizes = payload["size"]
    codes = [ordered_codes(payload["dimension"][dim]) for dim in ids]
    values = payload.get("value", {})
    rows = []
    for positions in itertools.product(*[range(size) for size in sizes]):
        flat = 0
        for pos, size in zip(positions, sizes):
            flat = flat * size + pos
        if isinstance(values, list):
            value = values[flat] if flat < len(values) else None
        else:
            value = values.get(str(flat), values.get(flat))
        row = {dim: codes[i][pos] for i, (dim, pos) in enumerate(zip(ids, positions))}
        row["value"] = value
        rows.append(row)
    return rows


def main(output):
    params = [
        ("lang", "en"),
        ("freq", "A"),
        ("unit", "YR"),
        ("sex", "T"),
        ("age", "Y65"),
    ]
    params.extend(("geo", geo) for geo in GEOS)
    params.extend(("time", year) for year in YEARS)
    response = requests.get(API, params=params, timeout=60)
    response.raise_for_status()
    rows = [
        row
        for row in flatten_jsonstat(response.json())
        if row.get("geo") in GEOS and row.get("time") in YEARS
    ]
    if len(rows) != len(GEOS) * len(YEARS):
        raise RuntimeError(f"Expected 8 values, got {len(rows)}")
    if any(row["value"] is None or not math.isfinite(float(row["value"])) for row in rows):
        raise RuntimeError("Eurostat returned a missing or non-finite value")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["geo", "year", "LE65", "source_url"])
        writer.writeheader()
        for row in sorted(rows, key=lambda item: (item["geo"], item["time"])):
            writer.writerow(
                {
                    "geo": row["geo"],
                    "year": row["time"],
                    "LE65": row["value"],
                    "source_url": response.url,
                }
            )
    print(output)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: fetch_eurostat_le65.py OUTPUT.csv")
    main(sys.argv[1])
