# Reproducible code for the INAR--rough-Heston manuscript

This directory contains the code needed for the numerical results in the
latest manuscript.  The original source tree is left unchanged.  The three
subdirectories separate the experiments by output type:

```text
code-revised/
├── INAR-tables/   INAR simulation, European/path-dependent tables, and checks
├── IV-slices/     multi-alpha implied-volatility slices and Table 2
└── IV-surface/    implied-volatility surface and ATM-skew figures
```

The C++ programs use only the C++17 standard library and POSIX threads.  The
Python scripts require Python 3 with `numpy`, `pandas`, and `matplotlib`;
`jupyter` is additionally needed for `assump_verify.ipynb`.

All commands below assume that they are run from the indicated subdirectory.
Use `g++` instead of `clang++` if preferred.

## 1. INAR tables

The core simulator is `INAR-tables/inar_cdq_all_options_v2.cpp`.  The table
driver includes that file and writes the tau-convergence, alpha-one European,
and alpha-one path-dependent CSV files.

```sh
cd code-revised/INAR-tables
clang++ -O3 -std=c++17 -pthread rerun_inar_tables.cpp -o rerun_inar_tables
./rerun_inar_tables \
  --paths 500000 \
  --output-dir ../../reruns/reproduced_inar_tables
```

The classical Heston Euler benchmark used for the path-dependent comparison
is run separately:

```sh
clang++ -O3 -std=c++17 -pthread classic_Heston_EM.cpp -o classic_Heston_EM
./classic_Heston_EM
```

The finite-`tau` transform and bounded-strip checks can be inspected or
executed with:

```sh
python3 -m jupyter notebook assump_verify.ipynb
```

The table driver uses the machine's available worker count.  The manuscript
rerun used 500,000 paths, seed `123456`, and 10 worker threads.  Monte Carlo
values should therefore be compared statistically across machines rather than
bit-for-bit.  The runtime columns in the manuscript are hardware-dependent
and should not be expected to match exactly.  The Euler column of the
path-dependent comparison comes from the separate `classic_Heston_EM.cpp`
run; the table driver's path-dependent CSV contains the INAR estimates and
confidence intervals.

## 2. Implied-volatility slices

`IV-slices/smile_overlay.cpp` includes `iv_slice_compare.cpp` and runs the
Integrated, INAR-CDQ-FFT, iVi, and HQE methods.  The following command
recreates the four-panel paper-style experiment with the manuscript settings:

```sh
cd code-revised/IV-slices
OUT=../../reruns/reproduced_iv_slices
mkdir -p "$OUT"

python3 run_smile_overlay_multi_alpha.py \
  --alphas 0.55 0.62 0.80 0.95 \
  --paths 16777216 \
  --threads 10 \
  --integrated-factors 3 \
  --output-dir "$OUT" \
  --binary "$OUT/smile_overlay" \
  --panel-output "$OUT/smile_overlay_panel.png"

python3 summarize_smile_results.py \
  "$OUT/smile_overlay_alpha_0_55_t_1_paths_16777216_paper_smile.csv" \
  "$OUT/smile_overlay_alpha_0_62_t_1_paths_16777216_paper_smile.csv" \
  "$OUT/smile_overlay_alpha_0_80_t_1_paths_16777216_paper_smile.csv" \
  "$OUT/smile_overlay_alpha_0_95_t_1_paths_16777216_paper_smile.csv" \
  --output "$OUT/summary.csv"
```

The runner also supports `--short` for the one-month slice and
`--inar-richardson` for the Richardson variant.  The plot script can be
called directly on any compatible set of CSV files.  The reported slice
errors are reproducible with the same seeds, path count, and thread count;
wall-clock times remain machine-dependent.

## 3. Implied-volatility surface and ATM skew

The surface generator reuses the INAR core through the relative include
`../INAR-tables/inar_cdq_all_options_v2.cpp`.  The default pipeline uses
1,000,000 paths, `tau=320`, 8 threads, seed `20260811`, and both rough and
classical models:

```sh
cd code-revised/IV-surface
./run_iv_surface_pipeline.sh
```

To write the output to a separate directory, override `IV_OUTDIR`:

```sh
IV_OUTDIR=../../reruns/reproduced_iv_surface \
IV_COMPILER=clang++ \
IV_PYTHON=python3 \
./run_iv_surface_pipeline.sh
```

The pipeline produces `surface_exact.csv`, `ivs.png`, `skew.png`,
`atm_diagnostics.csv`, and `summary.json`.

## Reference values and output files

The paper-style IV slice driver contains the precomputed rough-Heston
reference prices used for the comparison curves.  The tau-convergence
benchmark is likewise recorded in the table driver.  These values are part of
the reproducible input for the reported comparisons; an independent Fourier
reference-pricer is not required to run the Monte Carlo experiments.

For archival reproduction, retain the generated CSV files together with the
compiler, Python environment, path count, thread count, and random seeds.  Do
not rely on compiled binaries or Python `__pycache__` files as release
artifacts.
