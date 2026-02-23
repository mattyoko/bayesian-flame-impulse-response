function mcmc = runMCMC(logPostFunc, b0, propStruct, opts)
% Run a Metropolis-Hastings random-walk sampler without adaptation.
% Supports diagonal or full-covariance Gaussian proposals and returns both
% full chains and post-burn-in thinned samples.
%
% Inputs:
%   logPostFunc - Function handle: logPostFunc(b) -> scalar log posterior.
%   b0          - Initial parameter vector [Nb x 1].
%   propStruct  - Proposal definition:
%                 * vector of standard deviations, or
%                 * struct with .std and/or .cov (covariance takes precedence).
%   opts        - Optional struct with fields:
%                 nIter (default 20000), burnFrac (default 0.25),
%                 thin (default 1).
%
% Outputs:
%   mcmc - Struct containing chain history, acceptance diagnostics,
%          post-burn-in samples, and proposal metadata.

if nargin < 4, opts = struct(); end
if ~isfield(opts,'nIter'), opts.nIter = 20000; end
if ~isfield(opts,'burnFrac'), opts.burnFrac = 0.25; end
if ~isfield(opts,'thin'), opts.thin = 1; end

nIter = opts.nIter;
Nb = numel(b0);
chain = zeros(Nb,nIter);
logpost = -inf(1,nIter);
accepted = false(1,nIter);
windowAR = zeros(1,ceil(nIter/100));

% Determine proposal type
useFull = false;
if isstruct(propStruct)
    if isfield(propStruct,'cov') && ~isempty(propStruct.cov)
        propCov = propStruct.cov;
        useFull = true;
    elseif isfield(propStruct,'std')
        propStd = propStruct.std(:);
    else
        error('propStruct must contain either .cov or .std');
    end
else
    propStd = propStruct(:);
end

if useFull
    [L,p] = chol(propCov,'lower');
    if p>0, error('Proposal covariance not SPD'); end
else
    L = [];
end

% Init
chain(:,1) = b0(:);
logpost(1) = logPostFunc(chain(:,1));

% MH loop
for i = 2:nIter
    if useFull
        step = L*randn(Nb,1);
    else
        step = randn(Nb,1).*propStd;
    end
    prop = chain(:,i-1) + step;

    lp_curr = logpost(i-1);
    lp_prop = logPostFunc(prop);
    if ~isfinite(lp_prop)
        lp_prop = -inf;
    end

    alpha = min(1, exp(lp_prop - lp_curr));
    if rand < alpha
        chain(:,i) = prop;
        logpost(i) = lp_prop;
        accepted(i) = true;
    else
        chain(:,i) = chain(:,i-1);
        logpost(i) = lp_curr;
    end

    if mod(i,fix(nIter/10)) == 0
        fprintf('MCMC iteration %d of %d\n',i,nIter);
    end
end

% Window acceptance rates (for diagnostics only)
for w=1:numel(windowAR)
    lo = (w-1)*100+1; hi = min(w*100,nIter);
    windowAR(w) = mean(accepted(lo:hi));
end

% Outputs
burnIn = floor(opts.burnFrac*nIter);
thin = opts.thin;
idx = burnIn+1:thin:nIter;
postSamples = chain(:,idx);
postLogPost = logpost(idx);
acceptRate = mean(accepted(burnIn+1:end));

mcmc = struct();
mcmc.chain = chain;
mcmc.logpost = logpost;
mcmc.accepted = accepted;
mcmc.acceptRate = acceptRate;
mcmc.burnIn = burnIn;
mcmc.thin = thin;
mcmc.idx = idx;
mcmc.postSamples = postSamples;
mcmc.postLogPost = postLogPost;
if useFull
    mcmc.propCov = propCov;
else
    mcmc.propStd = propStd;
end
mcmc.windowAR = windowAR;
end