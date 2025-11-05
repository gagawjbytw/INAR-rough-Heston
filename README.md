# INAR-rough-Heston

Monte Carlo simulation implementation for pricing European/Asian/Lookback/Barrier options under the INAR(∞) framework, as described in our paper "Rough Heston Models as the scaling limit of bivariate heavy-tailed nearly unstable cumulative INAR(∞) processes".

## Overview

This repository contains the implementation of a Monte Carlo simulation framework for pricing various types of options using the INAR(∞) approach. The code demonstrates the connection between nearly unstable heavy-tailed INAR(∞) processes and the rough Heston model through numerical experiments.

## Requirements

- C++ compiler with C++17 support
- Standard library support
- CMake (optional, for build automation)

## Directory Structure

```
INAR-rough-Heston/
├── src/              # Source code files
├── scripts/          # Build and run scripts
└── README.md         # Documentation
```

## Building and Running

### Windows
```bash
# Using the provided batch script
scripts/rough_INAR.bat
```

### Linux/macOS
```bash
# Using the provided shell script
chmod +x scripts/build_and_run.sh
./scripts/rough_INAR.sh
```

## Features

The program calculates various option prices using Monte Carlo simulation, including:
- European call and put options
- Asian call and put options
- Lookback call and put options
- Barrier options (up-and-in call, down-and-out put)

Default parameters are set to match those used in the paper:
- Initial stock price (S₀) = 100.0
- Multiple strike prices around S₀
- Barrier levels at 110% and 90% of S₀

## Performance

The simulation is optimized for modern processors and includes features like:
- Efficient random number generation
- Path-wise calculation of multiple option types
- Comprehensive statistics including confidence intervals

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@article{cai2024scaling,
  title={Scaling limit of heavy tailed nearly unstable cumulative INAR ($$\backslash$infty $) processes and rough fractional diffusions},
  author={Wang, Yingli and Cai, Chunhao and He, Ping and Wang, QingHua},
  journal={arXiv preprint arXiv:2403.11773},
  year={2024}
}

@article{wang2025rough,
  title={Rough Heston Models as the scaling limit of bivariate heavy-tailed nearly unstable cumulative INAR($\infty$) processes},
  author={Wang, Yingli and Cui, Zhenyu and Zhu, Lingjiong},
  journal={arXiv preprint arXiv:2503.18259},
  year={2025}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
