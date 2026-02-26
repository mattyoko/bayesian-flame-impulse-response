# Overview

Bayesian inference tools for estimating the flame impulse response from input-output time-series data, with model-order comparison and uncertainty quantification. We restrict the search space to the family of distributed Gaussian time delays. For a given number of pulses, we use Bayesian parameter inference to find the most probable model parameters, given in the input-output data. We repeat this for several candidate models and use Bayesian model ranking to identify the most likely candidate model, given the data.

## Repository layout

- `src/` core inference and signal-processing functions
- `utils/` plotting and system-identification helpers
- `data/` example LES data from Eder et al. 2023 [1].
- `generateFigures.m` reproduce figures from Yoko & Polifke, "Bayesian inference of flame impulse responses".

[1] Alexander J. Eder, Camilo F. Silva, Matthias Haeringer, Johannes Kuhlmann, and Wolfgang Polifke. Incompressible versus compressible large eddy simulation for the identification of premixed flame dynamics. International Journal of Spray and Combustion Dynamics, 2023

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
Lref = 50e-3;       % Reference length, e.g. flame length
Uref = 11.3;        % Reference velocity, e.g. bulk flow speed
T_c = Lref / Uref;  % Reference time scale, e.g. convective time scale

% Load LES data
u = load(pathToVelocityData)
q = load(pathToHRRData)
t = load(pathToTimeData)

% Run model ranking, comparing orders 1 to 5
[h, posterior, modelRanking] = inferImpulseResponse(u, q, t, T_c, ...
    modelOrders=1:5, config=config);
```

## Reproducing the figures

This repo will evolve over time. To reproduce the figures in Yoko & Polifke 2026, make sure you check out the correct version:

```git
git clone https://github.com/mattyoko/bayesian-flame-impulse-response
cd <repo>
git checkout paper-v1.0
```
Then run the MATLAB script generateFigures.m

## Citation

If you use this code in academic work, please cite the accompanying manuscript:

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
