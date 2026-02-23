function [h, posterior, modelRanking, mcmc] = inferImpulseResponse(u, q, t, T_c, varargin)
% Infer an impulse response model from input-output time series data.
% Fits each candidate model order, ranks candidates with evidence-based
% metrics, and returns the best model's posterior and impulse response.
%
% Inputs:
%   u          - Input signal u(t) as a vector.
%   q          - Output signal q(t) as a vector.
%   t          - Time vector for u and q in seconds.
%   T_c        - Reference time scale used for non-dimensionalization.
%   varargin   - Optional name-value pairs:
%                * modelOrders: vector of candidate model orders.
%                * config: configuration struct from loadDefaultConfig.
%
% Outputs:
%   h            - Impulse response struct for the selected model.
%   posterior    - Posterior struct with MAP parameters and covariances.
%   modelRanking - Struct with per-model logML, logBFL, and logOF values.
%   mcmc         - MCMC output struct when enabled, otherwise empty.

% Parse optional inputs
p = inputParser;
addParameter(p, 'modelOrders', 1:5);
defaultConfig = loadDefaultConfig();
addParameter(p, 'config', defaultConfig);
parse(p, varargin{:}); 
modelOrders = p.Results.modelOrders;
config      = p.Results.config;
mcmc        = [];

% Initialize storage
names       = cell(size(modelOrders));
bMAP        = cell(size(modelOrders));
CbMAP       = cell(size(modelOrders));
hMAP        = cell(size(modelOrders));
CeMAP       = cell(size(modelOrders));
logML       = nan(size(modelOrders));
logBFL      = nan(size(modelOrders));
logOF       = nan(size(modelOrders));

Nmod = numel(modelOrders);

% Throw an error if the user has specified T_h but requested multiple model orders, 
% since this will bias model ranking
if ~isempty(config.model.T_h) && Nmod > 1
    error(['Model ranking across candidate orders is not meaningful for fixed T_h (see §4.4 of Yoko & Polifke 2026). ' ...
           'Either (i) set config.model.T_h = [], or (ii) only specify one model order in modelOrders.']);
end

% For each candidate model order, estimate the posterior:
col = getColours(1,Nmod);
for n = 1:Nmod
    order = modelOrders(n);
    % Generate prior
    [bp,Cb,names{n},t_max] = generatePrior(order,T_c,config);

    T_h = t_max * T_c; % Longest required impulse response signal
    % Prepare the signals:
    fs = 1/mean(diff(t));
    signals = prepareSignals(u, q, fs, T_h, config);
    % Call estimate posterior function, which minimizes the negative log posterior and then expands around the MAP estimate
    fprintf('Estimating posterior for %d delays: \n',order)
    [bMAP{n},CbMAP{n},hMAP{n},CeMAP{n},logML(n),logOF(n),logBFL(n)] = estimatePosterior(signals, bp, Cb, T_c, config);

    if config.diagnostics.plotCandidateHRRFit
        [~,~,~,~,~,p] = calculateCost(signals,CeMAP{n},bMAP{n},bp,Cb,T_c, config);
        figure(9001)
        if n == 1
            plot(signals.time_domain.coarse.t,signals.time_domain.coarse.q,'color',getColours(2),'lineWidth',2);
        end
        hold on;
        Nh = length(signals.time_domain.coarse.t_h);
        errorpatch(signals.time_domain.coarse.t(Nh:end),p,2*sqrt(CeMAP{n}),2*sqrt(CeMAP{n}),'color',col(n,:));
    end

    if config.diagnostics.plotCandidateImpulse
        figure(9002); hold on
        errorpatch(hMAP{n}.time,hMAP{n}.val,2*sqrt(hMAP{n}.var),2*sqrt(hMAP{n}.var),'color',col(n,:));
    end
end

modelRanking.logML = logML;
modelRanking.logBFL = logBFL;
modelRanking.logOF = logOF;

% Normalize model ranking to the most likely model
[logMLmax,indMax] = max(logML);
logML = logML - logMLmax;
logBFL = logBFL - logMLmax;

if numel(modelOrders) > 1
    % Print results to the command window
    fprintf('\n\nModel ranking outcome:\n\n')
    rows = compose('Delays: %d',modelOrders);
    table = array2table([logML(:) , logBFL(:) , logOF(:)],'rowNames',rows,'VariableNames',{'logML','logBFL','logOF'});
    format shortE
    disp(table);
    format short
    
    fprintf('A model with %d delays is most likely\n',modelOrders(indMax));
end

% Package for return
h       = hMAP{indMax};
h.time  = h.time .* T_c; % Undo time nondimensionalization
bMAP    = bMAP{indMax};
CbMAP   = CbMAP{indMax};
[a,~,~,Ca] = mapParameterToPhysicalSpace(bMAP,T_c,CbMAP);
Ce      = CeMAP{indMax};

posterior.a = a;
posterior.Ca = Ca;
posterior.b = bMAP;
posterior.Cb = CbMAP;
posterior.Ce = Ce;

if config.inference.runMCMC
    % Run MCMC for the best model
    % ------------------------ MCMC Sampling -----------------------
    [bp, Cb, paramNames] = generatePrior(modelOrders(indMax),T_c,config);
    Ce = CeMAP{indMax};
    logPosterior = @(b) -calculateCost(signals,Ce,b,bp,Cb,T_c,config); % returns log posterior up to constant
    % Use MAP covariance scaled for proposal
    d = numel(bMAP);
    scale = (2.38^2 / d);
    propCov = CbMAP * scale;
    optsMCMC = struct('nIter',config.inference.MCMCiter,'burnFrac',config.inference.MCMCburninFrac,'thin',config.inference.MCMCthin);
    mcmc = runMCMC(logPosterior, bMAP, struct('cov',propCov), optsMCMC);
    fprintf('MCMC acceptance rate (post burn-in): %.2f%%\n', 100*mcmc.acceptRate);
     
    if config.inference.plotMCMC
        fm = figure(1005);
        set(fm,'Position',[0 0 490 420]);
        cornerHeatmap(mcmc.postSamples', bMAP, CbMAP, paramNames, 40)
    end
end
end


