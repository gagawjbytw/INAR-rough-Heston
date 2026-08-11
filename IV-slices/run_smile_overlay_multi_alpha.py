import argparse
import subprocess
import sys
from pathlib import Path


DEFAULT_ALPHAS = [0.55, 0.62, 0.80, 0.95]


def fmt_alpha(alpha: float) -> str:
    return f"{alpha:.2f}".replace(".", "_")


def main():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--alphas", nargs="*", type=float, default=DEFAULT_ALPHAS)
    parser.add_argument("--paths", type=int, default=1 << 24)
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--short", action="store_true")
    parser.add_argument("--no-paper-smile", dest="paper_smile", action="store_false")
    parser.set_defaults(paper_smile=True)
    parser.add_argument("--inar-richardson", action="store_true")
    parser.add_argument("--integrated-factors", type=int)
    parser.add_argument("--binary", default="/tmp/smile_overlay")
    parser.add_argument("--output-dir", default="figures/multi_alpha")
    parser.add_argument("--panel-output")
    args = parser.parse_args()
    if args.threads <= 0:
        parser.error("--threads must be positive")

    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)

    binary = Path(args.binary)
    compile_cmd = [
        "g++", "-std=c++17", "-O3", "-pthread",
        str(script_dir / "smile_overlay.cpp"), "-o", str(binary)
    ]
    subprocess.run(compile_cmd, check=True)

    csv_paths = []
    for alpha in args.alphas:
        maturity_tag = "t_1_12" if args.short else "t_1"
        style_tag = "_paper_smile" if args.paper_smile else ""
        suffix = "_inar_rich" if args.inar_richardson else ""
        csv_path = root / f"smile_overlay_alpha_{fmt_alpha(alpha)}_{maturity_tag}_paths_{args.paths}{style_tag}{suffix}.csv"
        cmd = [
            str(binary),
            "--alpha",
            f"{alpha:.2f}",
            "--paths",
            str(args.paths),
            "--threads",
            str(args.threads),
            "--output",
            str(csv_path),
        ]
        if args.paper_smile:
            cmd.append("--paper-smile")
        if args.integrated_factors is not None:
            cmd.extend(["--integrated-factors", str(args.integrated_factors)])
        if args.short:
            cmd.append("--t-1-12")
        else:
            cmd.append("--t-1")
        if args.inar_richardson:
            cmd.append("--inar-richardson")
        subprocess.run(cmd, check=True)
        csv_paths.append(csv_path)

    panel_output = Path(args.panel_output) if args.panel_output else root / (
        ("smile_overlay_panel_t_1_12" if args.short else "smile_overlay_panel_t_1")
        + ("_paper_smile" if args.paper_smile else "")
        + ("_inar_rich" if args.inar_richardson else "")
        + ".png"
    )

    plot_cmd = [
        sys.executable,
        str(script_dir / "plot_smile_overlay_panel.py"),
        *[str(path) for path in csv_paths],
        "--output",
        str(panel_output),
        "--no-title",
    ]
    subprocess.run(plot_cmd, check=True)

    print("Generated panel from:")
    for path in csv_paths:
        print(path)
    print(panel_output)


if __name__ == "__main__":
    main()
