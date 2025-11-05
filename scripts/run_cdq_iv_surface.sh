#!/usr/bin/env bash

# Build and run the implied-volatility surface generator (CDQ version)
# and capture its CSV stdout into a file.

set -euo pipefail

echo "--- Implied volatility surface script (CDQ) ---"

# Move to the project root (script assumed to live in scripts/)
cd "$(dirname "$0")/.."

# Define file and directory paths
SRC_PRIMARY="src/cdq_mc_iv.cpp"
SRC_FALLBACK="cdq_mc_iv.cpp"
#SRC_PRIMARY="src/cdq_mc_iv.cpp"
#SRC_FALLBACK="cdq_mc_iv.cpp"
BUILD_DIR="build"
RESULTS_DIR="results"
EXECUTABLE_NAME="iv_surface_generator_cdq"
EXECUTABLE_PATH="$BUILD_DIR/$EXECUTABLE_NAME"
OUTPUT_CSV="$RESULTS_DIR/iv_surface_cdq.csv"

# Choose the source file
if [ -f "$SRC_PRIMARY" ]; then
  SRC_FILE="$SRC_PRIMARY"
elif [ -f "$SRC_FALLBACK" ]; then
  SRC_FILE="$SRC_FALLBACK"
else
  echo "Error: source file not found. Please store the C++ code as '$SRC_PRIMARY' or '$SRC_FALLBACK'."
  exit 1
fi

# Ensure g++ is available
if ! command -v g++ >/dev/null 2>&1; then
  echo "Error: g++ compiler not detected. Install it first (e.g., sudo apt-get install g++)."
  exit 1
fi

# Create build/results directories if missing
echo "Ensuring directories exist: '$BUILD_DIR' and '$RESULTS_DIR'..."
mkdir -p "$BUILD_DIR" "$RESULTS_DIR"

# Compile the C++ code
echo "Compiling '$SRC_FILE'..."
g++ -std=c++17 -O3 -march=native -Wall -Wextra -pthread "$SRC_FILE" -o "$EXECUTABLE_PATH"
echo "Compilation succeeded. Executable: $EXECUTABLE_PATH"

# Run the executable and redirect stdout to the CSV file
# Note: the program prints progress to stderr; stdout is pure CSV (with header)
echo "Running the model to generate surface data... this may take a while."
echo "Progress logs will appear below; final data goes to '$OUTPUT_CSV'"
echo "--------------------------------------------------------"

"$EXECUTABLE_PATH" > "$OUTPUT_CSV"

echo "--------------------------------------------------------"
echo "Done! Implied-volatility surface data saved to: $OUTPUT_CSV"
echo "--- Script finished ---"

# If you only want strict CSV (e.g., stripping any non-header lines),
# enable the block below. It backs up the raw output to *.raw before filtering:
# RAW="$OUTPUT_CSV.raw"
# mv "$OUTPUT_CSV" "$RAW"
# awk 'NR==1 && $0 !~ /^Maturity,/{next} {print $0}' "$RAW" > "$OUTPUT_CSV"
# echo "Filtered non-CSV lines; raw output stored as: $RAW"
