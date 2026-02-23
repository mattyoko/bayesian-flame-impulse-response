function config = loadDefaultConfig()
% Build the default configuration struct for inference and diagnostics.
% Centralizes preprocessing, optimizer, model, inference, and plotting
% settings used across the codebase.
%
% Inputs:
%   None.
%
% Outputs:
%   config - Struct containing default settings for all pipeline stages.

config = struct();

% Preprocessing settings
config.preproc.DSlimit = [];

% Prior settings    n    gam   beta
config.prior.mu  = [0.0, 0.0, -1.8];
config.prior.sig = [1.0, 0.5,  0.5];

% Optimizer settings
config.optimizer.restartScaling = 10;     % numRestarts = restartScaling * numParams
config.optimizer.iterationScaling = 100;  % maxIter = iterationScaling * numParams
config.optimizer.optimalityTolerance = 1e-5;
config.optimizer.stepSizeTolerance = 1e-10;
config.optimizer.functionChangeTolerance = 1e-10;
config.optimizer.finalHessianMethod = 'FD'; % 'GN' | 'BFGS' | 'FD'
config.optimizer.solnSelection = 'logProb'; % 'logML' | 'logProb'
config.optimizer.useParallel = true;

% Model settings
config.model.T_h = [];   % Total length of the impulse response. Blank calculates 95% bound from prior and model size
config.model.LFL = [];   % Low frequency limit ([] if free)
config.model.LFL_sigma = 0.01; % Uncertainty in the low frequency limit

% Inference settings
config.inference.dataVariance = [];     % Empty to infer noise by maximizing marginal likelihood
% -- MCMC settings
config.inference.runMCMC = false;
config.inference.plotMCMC = false;
config.inference.MCMCiter = 100000;     % Number of MCMC iterations
config.inference.MCMCburninFrac = 0.25; % Fraction of initial iterations to discard for burn-in
config.inference.MCMCthin = 5;          % Only accept every few samples to reduce chain correlation

% Diagnostics
config.diagnostics.plotCandidateImpulse = false; % Plot impulse response for each candidate model on same axes
config.diagnostics.plotCandidateHRRFit = false;  % Plot data and model HRR on same axes
end

