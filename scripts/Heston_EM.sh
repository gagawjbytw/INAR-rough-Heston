#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Get the project root directory (parent of scripts directory)
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

echo "Compiling Heston Model Program..."

# Create bin directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/bin"

# Change to project root directory
cd "$PROJECT_ROOT"

# Compile the C++ program
g++ -std=c++17 -O3 -pthread src/classic_Heston_EM.cpp -o bin/heston_model

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful!"
    echo "Running the program..."
    echo
    ./bin/heston_model
else
    echo "Compilation failed!"
    exit 1
fi

echo
echo "Program execution completed." 