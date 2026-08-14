#!/usr/bin/env python3
"""Sampled transform-strip diagnostics for the INAR--rough-Heston paper.

The script checks the parameter sets used in the numerical section on

    eta in {1.5, 2.0, 2.5},  xi in [-20, 20].

It reports three diagnostics:

1. the largest sampled modulus of the continuous fractional Riccati solution;
2. the largest coarse/fine time-grid discrepancy on the sampled rectangle;
3. the largest sampled modulus of the exact finite-tau log transform.

These are numerical diagnostics on compact grids.  They do not prove
analyticity or the unbounded-strip assumptions in the manuscript.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Parameters:
    label: str
    alpha: float
    gamma: float
    rho: float
    nu: float
    theta: float
    v0: float
    taus: tuple[int, ...]
    maturity: float = 1.0


def paper_parameter_sets() -> list[Parameters]:
    first = Parameters(
        label="rough-regime",
        alpha=0.62,
        gamma=0.1,
        rho=-0.681,
        nu=0.331,
        theta=0.3156,
        v0=0.0392,
        taus=(40, 80, 160, 320),
    )
    slice_common = dict(
        gamma=0.3,
        rho=-0.7,
        nu=1.0,
        theta=0.02 / 0.3,
        v0=0.02,
        taus=(160,),
    )
    slices = [
        Parameters(label="IV-slice", alpha=alpha, **slice_common)
        for alpha in (0.55, 0.62, 0.80, 0.95)
    ]
    return [first, *slices]


def riccati_rhs(h: complex, z: complex, p: Parameters) -> complex:
    return (
        0.5 * (z * z - z)
        + p.gamma * (p.rho * p.nu * z - 1.0) * h
        + 0.5 * (p.gamma * p.nu) ** 2 * h * h
    )


def solve_fractional_riccati_pc(
    z: complex,
    p: Parameters,
    steps: int,
    blowup_threshold: float = 1.0e10,
) -> np.ndarray:
    """Diethelm-type fractional Adams predictor-corrector discretization."""
    alpha = p.alpha
    dt = p.maturity / steps
    h = np.zeros(steps + 1, dtype=np.complex128)
    fvals = np.zeros(steps + 1, dtype=np.complex128)
    fvals[0] = riccati_rhs(0.0j, z, p)

    lags = np.arange(steps + 2, dtype=float)
    predictor_lag = (lags[1:] ** alpha) - (lags[:-1] ** alpha)
    corrector_lag = (
        lags[2:] ** (alpha + 1.0)
        + lags[:-2] ** (alpha + 1.0)
        - 2.0 * lags[1:-1] ** (alpha + 1.0)
    )
    predictor_scale = dt**alpha / math.gamma(alpha + 1.0)
    corrector_scale = dt**alpha / math.gamma(alpha + 2.0)

    for n in range(steps):
        predictor = predictor_scale * np.dot(
            fvals[: n + 1], predictor_lag[n::-1]
        )
        f_predictor = riccati_rhs(predictor, z, p)

        # The j=0 corrector weight is distinct from the lag weights.
        a0 = n ** (alpha + 1.0) - (n - alpha) * (n + 1) ** alpha
        history = a0 * fvals[0]
        if n >= 1:
            history += np.dot(fvals[1 : n + 1], corrector_lag[n - 1 :: -1])
        h[n + 1] = corrector_scale * (history + f_predictor)
        fvals[n + 1] = riccati_rhs(h[n + 1], z, p)

        if not np.isfinite(h[n + 1]) or abs(h[n + 1]) > blowup_threshold:
            raise FloatingPointError(
                f"continuous Riccati instability at z={z}, step={n + 1}/{steps}"
            )
    return h


def beta_from_rho(rho: float) -> float:
    coefficient = 1.0 - 2.0 * rho * rho
    discriminant = 4.0 - 4.0 * coefficient * coefficient
    roots = (
        (2.0 + math.sqrt(discriminant)) / (2.0 * coefficient),
        (2.0 - math.sqrt(discriminant)) / (2.0 * coefficient),
    )
    admissible = [root for root in roots if root > 1.0]
    if not admissible:
        raise ValueError(f"rho={rho} gives no beta>1")
    return admissible[0]


def base_kernel(alpha: float, horizon: int) -> np.ndarray:
    values = np.zeros(horizon + 1, dtype=float)
    inverse_gamma = 1.0 / math.gamma(1.0 - alpha)
    values[1] = 1.0 - inverse_gamma
    if horizon >= 2:
        n = np.arange(2, horizon + 1, dtype=float)
        values[2:] = inverse_gamma * ((n - 1.0) ** (-alpha) - n ** (-alpha))
    return values


def exact_discrete_log_transform(
    z: complex,
    p: Parameters,
    tau: int,
    blowup_threshold: float = 1.0e100,
) -> tuple[complex, float]:
    """Evaluate the exact finite-tau recursion and return log Phi and max |G|."""
    horizon = int(math.floor(tau * p.maturity))
    beta = beta_from_rho(p.rho)
    mu = (
        p.theta
        * (1.0 + beta * beta)
        / (p.gamma * p.nu * p.nu * (1.0 + beta) ** 2)
    )
    xi0 = p.v0 / p.theta
    a_tau = 1.0 - p.gamma * tau ** (-p.alpha)
    epsilon = 1.0 - a_tau
    mu_tau = mu * tau ** (p.alpha - 1.0)
    c_tau = math.sqrt(p.theta * epsilon / (2.0 * mu * tau**p.alpha))
    one_minus_exp = -math.expm1(-c_tau)
    d_tau = -math.log1p(-(one_minus_exp * one_minus_exp))

    phi = base_kernel(p.alpha, horizon)
    phi_tau = a_tau * phi
    q_plus = phi_tau / (1.0 + beta)
    q_minus = beta * phi_tau / (1.0 + beta)

    prefix = np.cumsum(phi_tau)
    baseline = np.zeros(horizon + 1, dtype=float)
    for n in range(1, horizon + 1):
        previous_sum = prefix[n - 1]
        baseline[n] = mu_tau + xi0 * mu_tau * (
            (1.0 - previous_sum) / epsilon - previous_sum
        )

    theta_plus = z * (c_tau - d_tau)
    theta_minus = -z * c_tau
    g = np.zeros(horizon + 1, dtype=np.complex128)
    max_g = 0.0
    for n in range(horizon + 1):
        plus_history = 0.0j
        minus_history = 0.0j
        if n >= 1:
            reversed_history = g[n - 1 :: -1]
            plus_history = np.dot(q_plus[1 : n + 1], reversed_history)
            minus_history = np.dot(q_minus[1 : n + 1], reversed_history)
        with np.errstate(over="raise", invalid="raise"):
            g[n] = (
                np.exp(theta_plus + plus_history)
                + np.exp(theta_minus + minus_history)
                - 2.0
            )
        max_g = max(max_g, float(abs(g[n])))
        if not np.isfinite(g[n]) or max_g > blowup_threshold:
            raise FloatingPointError(
                f"discrete transform instability at z={z}, tau={tau}, n={n}"
            )

    log_phi = np.dot(baseline[1:], g[horizon - 1 :: -1])
    if not np.isfinite(log_phi):
        raise FloatingPointError(f"nonfinite log Phi at z={z}, tau={tau}")
    return complex(log_phi), max_g


def check_parameter_set(
    p: Parameters,
    eta_values: np.ndarray,
    xi_values: np.ndarray,
    coarse_steps: int,
    fine_steps: int,
) -> dict[str, object]:
    if fine_steps % coarse_steps != 0:
        raise ValueError("fine_steps must be an integer multiple of coarse_steps")
    refinement_stride = fine_steps // coarse_steps

    rectangle_h_max = 0.0
    refinement_error = 0.0
    worst_h_z = 0.0j
    for eta in eta_values:
        for xi in xi_values:
            z = eta + 1j * xi
            coarse = solve_fractional_riccati_pc(z, p, coarse_steps)
            fine = solve_fractional_riccati_pc(z, p, fine_steps)
            local_max = float(np.max(np.abs(fine)))
            local_error = float(
                np.max(np.abs(fine[::refinement_stride] - coarse))
            )
            if local_max > rectangle_h_max:
                rectangle_h_max = local_max
                worst_h_z = z
            refinement_error = max(refinement_error, local_error)

    log_transform_max = 0.0
    recursion_max = 0.0
    worst_discrete = ""
    for tau in p.taus:
        for eta in eta_values:
            for xi in xi_values:
                z = eta + 1j * xi
                log_phi, max_g = exact_discrete_log_transform(z, p, tau)
                if abs(log_phi) > log_transform_max:
                    log_transform_max = float(abs(log_phi))
                    worst_discrete = f"tau={tau}, z={z}"
                recursion_max = max(recursion_max, max_g)

    return {
        "label": p.label,
        "alpha": p.alpha,
        "taus": ",".join(str(value) for value in p.taus),
        "h_rectangle_max": rectangle_h_max,
        "h_worst_z": str(worst_h_z),
        "h_refinement_sup": refinement_error,
        "log_transform_rectangle_max": log_transform_max,
        "g_rectangle_max": recursion_max,
        "discrete_worst": worst_discrete,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--R", type=float, default=20.0)
    parser.add_argument("--xi-points", type=int, default=41)
    parser.add_argument("--coarse-steps", type=int, default=1000)
    parser.add_argument("--fine-steps", type=int, default=2000)
    parser.add_argument(
        "--eta-values", type=float, nargs="+", default=(1.5, 2.0, 2.5)
    )
    parser.add_argument(
        "--output", type=Path, default=Path("transform_strip_diagnostics.csv")
    )
    args = parser.parse_args()

    if args.xi_points < 3 or args.xi_points % 2 == 0:
        raise ValueError("xi-points must be an odd integer at least 3")
    eta_values = np.asarray(args.eta_values, dtype=float)
    if not np.any(np.isclose(eta_values, 2.0)):
        raise ValueError("eta-values must include 2.0")
    xi_values = np.linspace(-args.R, args.R, args.xi_points)

    rows: list[dict[str, object]] = []
    for parameters in paper_parameter_sets():
        print(f"checking {parameters.label}, alpha={parameters.alpha:.2f}", flush=True)
        row = check_parameter_set(
            parameters,
            eta_values,
            xi_values,
            args.coarse_steps,
            args.fine_steps,
        )
        rows.append(row)
        print(
            "  "
            f"max|h|={row['h_rectangle_max']:.8g}, "
            f"refinement={row['h_refinement_sup']:.3g}, "
            f"max|log Phi^tau|={row['log_transform_rectangle_max']:.8g}",
            flush=True,
        )

    write_csv(args.output, rows)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
