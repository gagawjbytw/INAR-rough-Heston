#!/usr/bin/env bash
# Build and run the INAR CDQ all-options pricer and generate CSV outputs (macOS compatible)

set -euo pipefail

echo "--- INAR all-options pricing script (CDQ, macOS compatible) ---"

# Move to the project root (script lives inside scripts/)
cd "$(dirname "$0")/.."

# Paths
SRC_PRIMARY="src/inar_cdq_all_options.cpp"
SRC_FALLBACK="inar_cdq_all_options.cpp"
BUILD_DIR="build"
RESULTS_DIR="results"
EXECUTABLE_NAME="inar_all_options_cdq"
EXECUTABLE_PATH="$BUILD_DIR/$EXECUTABLE_NAME"

RAW_TXT="$RESULTS_DIR/inar_all_options_output.txt"
EU_CSV="$RESULTS_DIR/european_options.csv"
AS_CSV="$RESULTS_DIR/asian_options.csv"
LB_CSV="$RESULTS_DIR/lookback_options.csv"
BR_CSV="$RESULTS_DIR/barrier_options.csv"

# Pick the source file
if [ -f "$SRC_PRIMARY" ]; then
  SRC_FILE="$SRC_PRIMARY"
elif [ -f "$SRC_FALLBACK" ]; then
  SRC_FILE="$SRC_FALLBACK"
else
  echo "Error: source file not found. Please save it as '$SRC_PRIMARY' or '$SRC_FALLBACK'."
  exit 1
fi

# Ensure g++ is available
if ! command -v g++ >/dev/null 2>&1; then
  echo "Error: g++ was not detected. Install it first (e.g., xcode-select --install or brew install gcc)."
  exit 1
fi

# Directories
echo "Making sure '$BUILD_DIR' and '$RESULTS_DIR' exist..."
mkdir -p "$BUILD_DIR" "$RESULTS_DIR"

# Compile
echo "Compiling '$SRC_FILE' ..."
g++ -std=c++17 -O3 -march=native -Wall -Wextra -pthread "$SRC_FILE" -o "$EXECUTABLE_PATH"
echo "Build succeeded: $EXECUTABLE_PATH"

# Run and save raw output
echo "Running the program (progress/info prints to the terminal). Full output will be stored in: $RAW_TXT"
"$EXECUTABLE_PATH" | tee "$RAW_TXT" >/dev/null

# Parse stdout into CSV — strip [ ] and commas, then split on whitespace (BSD awk compatible)
echo "Parsing output into CSV ..."

awk -v euro_csv="$EU_CSV" \
    -v asian_csv="$AS_CSV" \
    -v look_csv="$LB_CSV" \
    -v barrier_csv="$BR_CSV" '
BEGIN{
  sec=""; # current section: EU/AS/LB/BR
  print "Strike,CallMean,CallCI_L,CallCI_U,PutMean,PutCI_L,PutCI_U" > euro_csv;
  print "Strike,CallMean,CallCI_L,CallCI_U,PutMean,PutCI_L,PutCI_U" > asian_csv;
  print "Strike,CallMean,CallCI_L,CallCI_U,PutMean,PutCI_L,PutCI_U" > look_csv;
  print "Strike,UpInCallMean,UpInCI_L,UpInCI_U,DownOutPutMean,DownOutCI_L,DownOutCI_U" > barrier_csv;
}
{
  # Detect section headers
  if ($0 ~ /^European Options:/) { sec="EU"; next }
  if ($0 ~ /^Asian Options:/)    { sec="AS"; next }
  if ($0 ~ /^Lookback Options:/) { sec="LB"; next }
  if ($0 ~ /^Barrier Options/)   { sec="BR"; next }

  # Skip separators, blank lines, and headers
  if ($0 ~ /^[-]+$/) next
  if ($0 ~ /^[[:space:]]*$/) next
  if ($0 ~ /^Strike[[:space:]]+Call/ || $0 ~ /^Strike[[:space:]]+Up-In-Call/) next

  # Only parse data rows that start with a number (Strike)
  if ($0 ~ /^[[:space:]]*[0-9]+[[:space:]]+/) {
    line=$0
    gsub(/\[/,"",line)
    gsub(/\]/,"",line)
    gsub(/,/,"",line)

    # Split on whitespace
    n=split(line,a)
    # Need at least 7 fields: Strike Mean L U Mean L U
    if (n>=7) {
      out=a[1]","a[2]","a[3]","a[4]","a[5]","a[6]","a[7]
      if (sec=="EU")       print out >> euro_csv
      else if (sec=="AS")  print out >> asian_csv
      else if (sec=="LB")  print out >> look_csv
      else if (sec=="BR")  print out >> barrier_csv
    }
  }
}
' "$RAW_TXT"

echo "CSV files generated:"
echo "  - $EU_CSV"
echo "  - $AS_CSV"
echo "  - $LB_CSV"
echo "  - $BR_CSV"
echo "--- Script completed ---"
