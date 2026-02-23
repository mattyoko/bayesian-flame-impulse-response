function [bp,Cp,names,t_max] = generatePrior(N,T_c,config)
% Construct the prior over model parameters for a chosen model order.
% Returns Gaussian priors in parameter space and an estimated maximum
% non-dimensional horizon used for signal preparation.
%
% Inputs:
%   N      - Number of pulse components (model order).
%   T_c    - Reference time scale used for non-dimensionalization.
%   config - Configuration struct (uses config.model.T_h when provided).
%
% Outputs:
%   bp     - Prior mean vector in parameter space [3N x 1].
%   Cp     - Prior covariance matrix in parameter space [3N x 3N].
%   names  - Parameter label vector for plotting/reporting.
%   t_max  - Maximum non-dimensional impulse-response horizon.

% Create parameter name vector
names = compose('$%s_%d$',["n","\gamma","\beta"]',1:N);
names = names(:);

% Extract per-parameter prior from config
mu = config.prior.mu(:);
sig= config.prior.sig(:);

% Create prior mean vector
bp = repmat(mu,1,N);
bp = bp(:);

% Create prior covariance matrix
sigp = repmat(sig,1,N);
sigp = sigp(:);
Cp = diag(sigp.^2);

if isempty(config.model.T_h)
    % Calculate the 99% confidence upper bound on the length of the impulse
    % response, given this prior
    n99 = 2.326; % Number of standard deviations to reach 99% confidence
    
    % Calculate Fenton-Wilkinson estimate for longest probable sum of
    % log-normally distributed delay gaps
    mutau  = bp(2);
    Vtau   = sigp(2)^2;
    % Expectation and variance of the sum:
    ES = N * exp(mutau + 0.5 * Vtau);
    VS = N * (exp(Vtau) - 1) * exp(2*mutau + Vtau);
    % Fenton-Wilkinson estimate:
    sigFW  = sqrt(log(1 + VS/ES^2));
    muFW   = log(ES) - 0.5*sigFW^2;
    S99    = exp(muFW + n99 * sigFW);
    
    % Add longest possible spread
    musig   = bp(3);
    Vsig    = sigp(3)^2;
    sig_max = exp(musig +n99*sqrt(Vsig));
    
    t_max = S99 + n99*sig_max;
else
    t_max = config.model.T_h / T_c;
end

end