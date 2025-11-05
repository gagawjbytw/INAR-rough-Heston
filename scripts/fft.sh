#!/bin/bash

set -e

echo "--- FFT options pricer: build and run script ---"

cd "$(dirname "$0")/.."

SRC_FILE="src/fft.cpp"
BUILD_DIR="build"
EXECUTABLE_NAME="fft_pricer"
EXECUTABLE_PATH="$BUILD_DIR/$EXECUTABLE_NAME"

if [ ! -f "$SRC_FILE" ]; then
    echo "Error: source file '$SRC_FILE' not found!"
    exit 1
fi

echo "Creating build directory '$BUILD_DIR'..."
mkdir -p "$BUILD_DIR"

echo "Compiling '$SRC_FILE'..."
g++ -std=c++17 -O3 -Wall -pthread "$SRC_FILE" -o "$EXECUTABLE_PATH"

if [ $? -eq 0 ]; then
    echo "Compilation succeeded!"
else
    echo "Compilation failed."
    exit 1
fi

echo "Running executable '$EXECUTABLE_PATH'..."
echo "----------------------------------------------"
"$EXECUTABLE_PATH"
echo "----------------------------------------------"
echo "Program finished."
echo "--- Script finished ---"