import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


METHOD_STYLE = {
    "Reference": {
        "label": "Reference",
        "color": "black",
        "linestyle": "-",
        "marker": "o",
        "linewidth": 2.0,
        "markersize": 4.0,
    },
    "Integrated-MultifactorEuler": {
        "label": "Integrated",
        "color": "#8c564b",
        "linestyle": "-",
        "marker": "s",
        "linewidth": 1.7,
        "markersize": 4.0,
    },
    "INAR-CDQ-FFT": {
        "label": "INAR",
        "color": "#d62728",
        "linestyle": "-",
        "marker": "x",
        "linewidth": 1.7,
        "markersize": 4.2,
    },
    "INAR-Richardson(p=1)": {
        "label": "INAR-Richardson",
        "color": "#ff7f0e",
        "linestyle": "--",
        "marker": "D",
        "linewidth": 1.7,
        "markersize": 3.8,
    },
    "iVi-InverseGaussian": {
        "label": "iVi",
        "color": "#1f77b4",
        "linestyle": "-",
        "marker": "o",
        "linewidth": 1.7,
        "markersize": 4.2,
    },
    "HQE-RoughHeston": {
        "label": "HQE",
        "color": "#2ca02c",
        "linestyle": "-",
        "marker": "^",
        "linewidth": 1.7,
        "markersize": 4.2,
    },
}

METHOD_ORDER = [
    "Reference",
    "Integrated-MultifactorEuler",
    "INAR-CDQ-FFT",
    "INAR-Richardson(p=1)",
    "iVi-InverseGaussian",
    "HQE-RoughHeston",
]


def read_rows(csv_path: Path):
    grouped = {}
    alpha = None
    maturity = None
    paths = None
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            method = row["method"]
            grouped.setdefault(method, {"x": [], "iv": []})
            grouped[method]["x"].append(float(row["log_moneyness"]))
            grouped[method]["iv"].append(float(row["iv"]))
            alpha = float(row["alpha"])
            maturity = float(row["T"])
            paths = int(row["paths"])

    for payload in grouped.values():
        pairs = sorted(zip(payload["x"], payload["iv"]))
        payload["x"] = [x for x, _ in pairs]
        payload["iv"] = [iv for _, iv in pairs]

    if alpha is None or maturity is None or paths is None:
        raise RuntimeError(f"No rows found in {csv_path}")
    return grouped, alpha, maturity, paths


def maturity_label(maturity: float) -> str:
    if abs(maturity - 1.0) < 1e-12:
        return "T = 1"
    if abs(maturity - (1.0 / 12.0)) < 1e-12:
        return "T = 1/12"
    return f"T = {maturity:g}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_paths", nargs="+")
    parser.add_argument("--output")
    parser.add_argument("--no-title", action="store_true")
    args = parser.parse_args()

    csv_paths = [Path(p) for p in args.csv_paths]
    if len(csv_paths) == 0:
        raise RuntimeError("Need at least one CSV.")
    if len(csv_paths) > 4:
        raise RuntimeError("At most 4 CSVs are supported for a 2x2 panel.")

    payloads = [read_rows(path) for path in csv_paths]
    payloads.sort(key=lambda item: item[1])

    _, _, maturity0, paths0 = payloads[0]
    for _, _, maturity, paths in payloads[1:]:
        if abs(maturity - maturity0) > 1e-12:
            raise RuntimeError("All panel CSVs must share the same maturity.")
        if paths != paths0:
            raise RuntimeError("All panel CSVs must share the same paths count.")

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 9,
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.8), sharex=True, sharey=True)
    axes = axes.flatten()

    for ax_idx, (grouped, alpha, _, _) in enumerate(payloads):
        ax = axes[ax_idx]
        for method in METHOD_ORDER:
            if method not in grouped:
                continue
            style = METHOD_STYLE[method]
            ax.plot(
                grouped[method]["x"],
                grouped[method]["iv"],
                label=style["label"],
                color=style["color"],
                linestyle=style["linestyle"],
                marker=style["marker"],
                linewidth=style["linewidth"],
                markersize=style["markersize"],
            )
        ax.set_title(rf"$\alpha={alpha:.2f}$")
        ax.grid(True, color="#d0d0d0", linewidth=0.6, alpha=0.7)

    for ax in axes[len(payloads):]:
        ax.axis("off")

    axes[0].set_ylabel("Implied volatility")
    axes[2].set_ylabel("Implied volatility")
    axes[2].set_xlabel("Log-moneyness")
    axes[3].set_xlabel("Log-moneyness")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        frameon=True,
        facecolor="white",
        edgecolor="black",
        framealpha=1.0,
        bbox_to_anchor=(0.5, 0.98),
    )

    if not args.no_title:
        title = (
            rf"Implied-volatility smile overlays, {maturity_label(maturity0)}, "
            + (rf"paths=$2^{{24}}$" if paths0 == (1 << 24) else f"paths={paths0:,}")
        )
        fig.suptitle(title, y=0.995)
        rect = [0, 0, 1, 0.95]
    else:
        rect = [0, 0, 1, 0.93]

    fig.tight_layout(rect=rect)

    if args.output:
        out_path = Path(args.output)
    else:
        stem = f"smile_overlay_panel_{maturity_label(maturity0).replace(' ', '').replace('=', '').replace('/', '_')}"
        out_path = csv_paths[0].with_name(stem + ".png")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Saved {out_path}")
    print(f"Saved {out_path.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
