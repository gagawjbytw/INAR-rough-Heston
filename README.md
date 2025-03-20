# INAR-rough-Heston

Monte Carlo simulation implementation for pricing European options under the INAR(∞) framework, as described in our paper "Rough Heston Models as the scaling limit of bivariate heavy-tailed nearly unstable cumulative INAR(∞) processes".

## Overview

This repository contains the implementation of a parallel Monte Carlo simulation framework for pricing European options using the INAR(∞) approach. The code demonstrates the connection between nearly unstable heavy-tailed INAR(∞) processes and the rough Heston model through numerical experiments.

## Requirements

- C++ compiler with C++11 support
- pthread library support
- CMake (optional, for build automation)

## Directory Structure

```
INAR-rough-Heston/
├── src/              # Source code files
├── scripts/          # Build and run scripts
├── examples/         # Example usage and test cases
└── doc/             # Documentation
```

## Building and Running

### Windows
```bash
# Using the provided batch script
scripts/build_and_run.bat

# Or manually
g++ src/monte_carlo_inar_multi_strikes.cpp -o bin/monte_carlo_inar_multi_strikes.exe -std=c++11 -pthread -O2
```

## Usage Example

The program calculates European call option prices for multiple strikes using parallel Monte Carlo simulation. Default parameters are set to match those used in the paper:

- Hurst parameter (α) = 0.62
- Scaling parameter (γ) = 0.1
- Initial stock price (S₀) = 100.0
- Number of time steps = 1000
- Number of Monte Carlo paths = 100,000

## Performance

On a modern multi-core processor, the simulation typically completes within a few minutes, with the exact time depending on the number of available CPU cores and the chosen parameters.

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@article{your-paper,
  title={Rough Heston Models as the scaling limit of bivariate heavy-tailed nearly unstable cumulative INAR(∞) processes},
  author={Cui, Zhenyu and Wang, Yingli},
  journal={},
  year={2024}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. 