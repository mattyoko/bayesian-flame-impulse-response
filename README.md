# ImpulseResponse

Bayesian inference tools for estimating an impulse response from input-output time-series data, with model-order comparison and uncertainty quantification.

## Overview

This repository contains MATLAB code to:
- infer impulse responses using a pulse-based parameterization,
- estimate posterior uncertainty via Laplace approximation,
- optionally sample posteriors using MCMC,
- compare candidate model orders with evidence-based scores,
- compute and visualize transfer-function quantities.

## Repository layout

- `src/` core inference and signal-processing functions
- `utils/` plotting and system-identification helpers
- `data/` example LES data from Eder et al. 2023.
- `generateFigures.m` reproduce figures from Yoko & Polifke, "Bayesian inference of flame impulse responses".

## Getting started

1. Open the repository in MATLAB.
2. Add `src/`, `data/` and `utils/` to the MATLAB path.
3. Load or prepare input/output signals `u(t)` and `q(t)` with time vector `t`.
4. Run inference via `inferImpulseResponse(...)`.

Minimal example:

```matlab
% Set the default config, adjust however you want
config = loadDefaultConfig();
% Set the convective timescale of the system
Lref = 50e-3;
Uref = 11.3;
T_c = Lref / Uref;

% Load LES data
u = load(pathToVelocityData)
q = load(pathToHRRData)
t = load(pathToTimeData)

% Run model ranking, comparing orders 1 to 5
[h, posterior, modelRanking] = inferImpulseResponse(u, q, t, T_c, ...
    modelOrders=1:5, config=config);
```

## Citation

If you use this code in academic work, please cite the accompanying manuscript (preferred):

- Matthew Yoko and Wolfgang Polifke, “Bayesian Inference of Flame Impulse Responses”, manuscript under review, 2026.

BibTeX:

```bibtex
@unpublished{YokoPolifke2026ImpulseResponse,
    title  = {Bayesian Inference of Flame Impulse Responses},
    author = {Yoko, Matthew and Polifke, Wolfgang},
    year   = {2026},
    note   = {Manuscript under review},
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
