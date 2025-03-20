# Technical Documentation

## Implementation Details

This document provides technical details about the implementation of the INAR(∞) Monte Carlo simulation for option pricing.

### Core Components

1. **Fractional Kernel Function (`f_alpha_one`)**
   - Implements the fractional kernel \[f_\alpha(t)\] as defined in the paper
   - Handles special cases for t ≤ 0 and t = 1
   - Uses the gamma function for general case calculations

2. **Path Simulation (`simulate_inar_path`)**
   - Generates a single price path using the INAR(∞) process
   - Key steps:
     - Computes scaling parameters a_T and μ_T
     - Initializes and maintains counting processes (X⁺, X⁻, N⁺, N⁻)
     - Computes intensities λ⁺ and λ⁻
     - Generates Poisson random variables for order flow
     - Computes final price through exponential transformation

3. **Parallel Monte Carlo Implementation**
   - Uses C++11 threading facilities
   - Divides simulations across available CPU cores
   - Thread-safe random number generation
   - Efficient memory management through vector pre-allocation

### Performance Considerations

1. **Memory Management**
   - Pre-allocated vectors for phi and phi_cum
   - Reuse of vectors for counting processes
   - Efficient handling of large arrays

2. **Computational Optimization**
   - Parallel execution using std::async
   - Thread-local random number generators
   - Vectorized operations where possible

3. **Numerical Stability**
   - Careful handling of exponential terms
   - Bounds checking for Poisson parameters
   - Accumulation of sums using appropriate precision

### Parameter Settings

The default parameters are chosen to match those in the paper:
- α = 0.62 (Hurst parameter)
- γ = 0.1 (Scaling parameter)
- β = 27.558 (Asymmetry parameter)
- μ = 26.859152 (Base intensity)
- θ = 0.3156 (Volatility scaling)
- ξ = 0.1242 (Intensity adjustment)

### Validation

The implementation has been validated against:
1. Results from the paper
2. Alternative numerical methods (Hybrid Schemes and Adams scheme)
3. Known limiting behaviors

## Usage Guidelines

1. **Compilation**
   ```bash
   g++ src/monte_carlo_inar_multi_strikes.cpp -o bin/monte_carlo_inar_multi_strikes.exe -std=c++11 -pthread -O2
   ```

2. **Parameter Adjustment**
   - Parameters can be modified in the main function
   - Example configuration file provided in examples/example_config.txt

3. **Output Format**
   - Strike prices and corresponding option prices
   - Execution time information
   - Precision set to 4 decimal places

## Future Improvements

1. **Potential Optimizations**
   - GPU acceleration for larger simulations
   - Improved memory management for very large T_steps
   - Advanced variance reduction techniques

2. **Additional Features**
   - Configuration file support
   - More option types
   - Real-time progress reporting
   - Parameter calibration tools 