#!/usr/bin/env python3
"""Plot the paper IV-surface and ATM-skew diagnostics from generated CSV data.

Example:
    /opt/anaconda3/bin/python3 code/plot_iv_surface.py \
        --input code/iv_surface_output/surface_exact.csv \
        --outdir code/iv_surface_output

The script writes ``ivs.png``, ``skew.png``, ``atm_diagnostics.csv``, and
``summary.json``.  It never edits the manuscript or copies figures into a
submission directory automatically.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--dpi", default=220, type=int)
    return parser.parse_args()


def validate(data: pd.DataFrame) -> None:
    required = {
        "model",
        "alpha",
        "maturity",
        "log_moneyness",
        "implied_vol",
        "spot_mean",
        "spot_se",
        "paths",
        "tau",
        "threads",
        "c_tau",
        "d_tau",
        "negative_lambda_clamps",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"missing CSV columns: {sorted(missing)}")
    if set(data["model"].unique()) != {"rough", "classical"}:
        raise ValueError("the comparison plot requires both rough and classical rows")
    if data["implied_vol"].isna().any() or not np.isfinite(data["implied_vol"]).all():
        bad = data.loc[~np.isfinite(data["implied_vol"])]
        raise ValueError(f"non-finite implied vols in {len(bad)} rows")
    for model, group in data.groupby("model"):
        expected = group["maturity"].nunique() * group["log_moneyness"].nunique()
        if len(group) != expected:
            raise ValueError(f"{model}: surface grid is incomplete or has duplicate rows")


def matrix(data: pd.DataFrame, model: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    subset = data.loc[data["model"] == model]
    pivot = subset.pivot(index="log_moneyness", columns="maturity", values="implied_vol")
    k = pivot.index.to_numpy(dtype=float)
    maturity = pivot.columns.to_numpy(dtype=float)
    return maturity, k, pivot.to_numpy(dtype=float)


def plot_surfaces(data: pd.DataFrame, destination: Path, dpi: int) -> None:
    rough_t, rough_k, rough_iv = matrix(data, "rough")
    classical_t, classical_k, classical_iv = matrix(data, "classical")
    if not np.allclose(rough_t, classical_t) or not np.allclose(rough_k, classical_k):
        raise ValueError("rough and classical grids differ")

    maturity_mesh, k_mesh = np.meshgrid(rough_t, rough_k)
    figure = plt.figure(figsize=(9.3, 8.8))
    axis = figure.add_subplot(111, projection="3d")
    axis.plot_surface(
        maturity_mesh,
        k_mesh,
        rough_iv,
        cmap="viridis",
        alpha=0.72,
        linewidth=0.25,
        edgecolor=(1, 1, 1, 0.22),
        antialiased=True,
    )
    axis.plot_surface(
        maturity_mesh,
        k_mesh,
        classical_iv,
        cmap="plasma",
        alpha=0.55,
        linewidth=0.25,
        edgecolor=(1, 1, 1, 0.18),
        antialiased=True,
    )
    axis.set_title("Implied Volatility Surface (Comparison)", pad=22)
    axis.set_xlabel("Maturity (T)", labelpad=10)
    axis.set_ylabel("Log-Moneyness (k)", labelpad=10)
    axis.set_zlabel("Implied Volatility (%)", labelpad=8)
    axis.zaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0, decimals=1))
    axis.view_init(elev=29, azim=-122)
    axis.legend(
        handles=[
            Patch(facecolor=plt.cm.viridis(0.55), edgecolor="black", label="Rough (alpha = 0.62)"),
            Patch(facecolor=plt.cm.plasma(0.55), edgecolor="black", label="Classical (alpha = 1)"),
        ],
        loc="upper right",
        bbox_to_anchor=(1.03, 0.94),
    )
    figure.tight_layout()
    figure.savefig(destination, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def atm_diagnostics(data: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    rough = data.loc[data["model"] == "rough"].copy()
    k_values = np.sort(rough["log_moneyness"].unique())
    zero_index = int(np.argmin(np.abs(k_values)))
    if not np.isclose(k_values[zero_index], 0.0, atol=1e-12):
        raise ValueError("log-moneyness grid has no ATM point")
    if zero_index == 0 or zero_index == len(k_values) - 1:
        raise ValueError("ATM point has no central-difference neighbours")
    k_minus = float(k_values[zero_index - 1])
    k_plus = float(k_values[zero_index + 1])
    if not np.isclose(k_plus, -k_minus, atol=1e-12):
        raise ValueError("ATM neighbouring points are not symmetric")

    rows: list[dict[str, float]] = []
    for maturity, group in rough.groupby("maturity", sort=True):
        iv_by_k = group.set_index("log_moneyness")["implied_vol"]
        atm_iv = float(iv_by_k.loc[0.0])
        skew = float((iv_by_k.loc[k_plus] - iv_by_k.loc[k_minus]) / (k_plus - k_minus))
        first = group.iloc[0]
        rows.append(
            {
                "maturity": float(maturity),
                "atm_iv": atm_iv,
                "atm_skew": skew,
                "abs_atm_skew": abs(skew),
                "spot_mean": float(first["spot_mean"]),
                "spot_se": float(first["spot_se"]),
            }
        )
    diagnostics = pd.DataFrame(rows).sort_values("maturity").reset_index(drop=True)

    x = np.log(diagnostics["maturity"].to_numpy())
    y = np.log(diagnostics["abs_atm_skew"].to_numpy())
    slope, intercept = np.polyfit(x, y, 1)
    fitted = intercept + slope * x
    residual_sum = float(np.sum((y - fitted) ** 2))
    total_sum = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - residual_sum / total_sum if total_sum > 0.0 else math.nan
    summary = {
        "fit_exponent_slope": float(slope),
        "fit_beta_minus_slope": float(-slope),
        "fit_prefactor": float(np.exp(intercept)),
        "fit_r_squared": r_squared,
        "central_difference_k_minus": k_minus,
        "central_difference_k_plus": k_plus,
        "atm_iv_T_1": float(diagnostics.iloc[-1]["atm_iv"]),
        "atm_skew_T_1_12": float(diagnostics.iloc[0]["atm_skew"]),
        "atm_skew_T_1": float(diagnostics.iloc[-1]["atm_skew"]),
    }
    diagnostics["power_law_fit"] = np.exp(intercept) * diagnostics["maturity"] ** slope
    return diagnostics, summary


def plot_skew(
    diagnostics: pd.DataFrame,
    summary: dict[str, float],
    destination: Path,
    dpi: int,
) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(14.1, 4.25))
    maturity = diagnostics["maturity"].to_numpy()
    skew = diagnostics["atm_skew"].to_numpy()
    atm_iv = diagnostics["atm_iv"].to_numpy()
    absolute_skew = diagnostics["abs_atm_skew"].to_numpy()

    axes[0].plot(maturity, skew, marker="o", linewidth=1.6)
    axes[0].set_title("ATM Skew vs Maturity")
    axes[0].set_xlabel("Maturity (years)")
    axes[0].set_ylabel("ATM skew (d sigma / d log K)")
    axes[0].grid(alpha=0.35)

    axes[1].plot(maturity, atm_iv, marker="o", linewidth=1.6)
    axes[1].set_title("ATM IV vs Maturity")
    axes[1].set_xlabel("Maturity (years)")
    axes[1].set_ylabel("ATM implied volatility")
    axes[1].grid(alpha=0.35)

    order = np.argsort(maturity)
    axes[2].loglog(maturity, absolute_skew, "o", label="Monte Carlo")
    axes[2].loglog(
        maturity[order],
        diagnostics["power_law_fit"].to_numpy()[order],
        linewidth=1.6,
        label=(
            f"fit: beta={summary['fit_beta_minus_slope']:.3f}, "
            f"R^2={summary['fit_r_squared']:.3f}"
        ),
    )
    axes[2].set_title("Power-law fit (log scale)")
    axes[2].set_xlabel("Maturity (years) [log]")
    axes[2].set_ylabel("|ATM skew| [log]")
    axes[2].grid(which="both", alpha=0.35)
    axes[2].legend(frameon=False, fontsize=9)

    for index, axis in enumerate(axes):
        axis.text(0.02, -0.18, f"({chr(ord('a') + index)})", transform=axis.transAxes)
    figure.tight_layout()
    figure.savefig(destination, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input)
    validate(data)

    diagnostics, summary = atm_diagnostics(data)
    plot_surfaces(data, args.outdir / "ivs.png", args.dpi)
    plot_skew(diagnostics, summary, args.outdir / "skew.png", args.dpi)
    diagnostics.to_csv(args.outdir / "atm_diagnostics.csv", index=False)

    first_rows = data.groupby("model", sort=True).first()
    input_sha256 = hashlib.sha256(args.input.read_bytes()).hexdigest()
    maturity_steps = (
        data.loc[data["model"] == "rough", ["maturity", "step"]]
        .drop_duplicates()
        .sort_values("maturity")
    )
    summary.update(
        {
            "input": str(args.input),
            "input_sha256": input_sha256,
            "paths": int(data["paths"].iloc[0]),
            "tau_per_year": int(data["tau"].iloc[0]),
            "threads": int(data["threads"].iloc[0]),
            "rough_seed": int(first_rows.loc["rough", "seed"]),
            "classical_seed": int(first_rows.loc["classical", "seed"]),
            "gamma": float(data["gamma"].iloc[0]),
            "theta": float(data["theta"].iloc[0]),
            "rho": float(data["rho"].iloc[0]),
            "nu": float(data["nu"].iloc[0]),
            "V0": float(data["V0"].iloc[0]),
            "S0": float(data["S0"].iloc[0]),
            "beta": float(data["beta"].iloc[0]),
            "mu": float(data["mu"].iloc[0]),
            "xi": float(data["xi"].iloc[0]),
            "maturity_steps": maturity_steps.to_dict(orient="records"),
            "log_moneyness_grid": sorted(float(value) for value in data["log_moneyness"].unique()),
            "rough_spot_mean_T_1": float(
                data.loc[(data["model"] == "rough") & np.isclose(data["maturity"], 1.0), "spot_mean"].iloc[0]
            ),
            "classical_spot_mean_T_1": float(
                data.loc[(data["model"] == "classical") & np.isclose(data["maturity"], 1.0), "spot_mean"].iloc[0]
            ),
            "rough_c_tau": float(first_rows.loc["rough", "c_tau"]),
            "rough_d_tau": float(first_rows.loc["rough", "d_tau"]),
            "classical_c_tau": float(first_rows.loc["classical", "c_tau"]),
            "classical_d_tau": float(first_rows.loc["classical", "d_tau"]),
            "rough_seconds": float(first_rows.loc["rough", "seconds"]),
            "classical_seconds": float(first_rows.loc["classical", "seconds"]),
            "negative_lambda_clamps": int(data["negative_lambda_clamps"].max()),
        }
    )
    for model in ("rough", "classical"):
        c_tau = float(first_rows.loc[model, "c_tau"])
        d_tau = float(first_rows.loc[model, "d_tau"])
        summary[f"{model}_rn_identity_error"] = float(
            np.exp(c_tau - d_tau) + np.exp(-c_tau) - 2.0
        )
    with (args.outdir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
