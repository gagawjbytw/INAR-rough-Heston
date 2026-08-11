#!/usr/bin/env python3
"""Summarize multi-alpha smile CSV files into one reproducible table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DIAGNOSTIC_FIELDS = [
    "discounted_spot_mean",
    "spot_se",
    "martingale_z_score",
    "terminal_integrated_variance_mean",
    "terminal_noise_bracket_mean",
    "integrated_bracket_gap_max",
    "integrated_driver_mean",
    "integrated_driver_second_moment",
    "integrated_driver_bracket_ratio",
    "integrated_variance_driver_second_moment",
    "integrated_variance_driver_bracket_ratio",
    "projection_hit_rate",
    "projection_gap_mean",
    "projection_gap_max",
]


def summarize(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    references = {
        row["log_moneyness"]: float(row["iv"])
        for row in rows
        if row["method"] == "Reference"
    }
    if not references:
        raise RuntimeError(f"No Reference rows in {path}")

    result: list[dict[str, str]] = []
    methods = sorted({row["method"] for row in rows if row["method"] != "Reference"})
    for method in methods:
        method_rows = [row for row in rows if row["method"] == method]
        if set(row["log_moneyness"] for row in method_rows) != set(references):
            raise RuntimeError(f"Grid mismatch for {method} in {path}")
        maximum_error = max(
            abs(float(row["iv"]) - references[row["log_moneyness"]])
            for row in method_rows
        )
        first = method_rows[0]
        summary_row = {
            "alpha": first["alpha"],
            "T": first["T"],
            "paths": first["paths"],
            "threads": first.get("threads", ""),
            "method": method,
            "steps": first["steps"],
            "factors": first["factors"],
            "max_abs_iv_error": f"{maximum_error:.9g}",
            "time_seconds": first["time_seconds"],
            "source_csv": str(path),
        }
        for field in DIAGNOSTIC_FIELDS:
            summary_row[field] = first.get(field, "")
        result.append(summary_row)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_paths", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    for path in args.csv_paths:
        rows.extend(summarize(path))
    rows.sort(key=lambda row: (float(row["alpha"]), row["method"]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "alpha", "T", "paths", "threads", "method", "steps", "factors",
        "max_abs_iv_error", "time_seconds", *DIAGNOSTIC_FIELDS, "source_csv",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
