function [J,dJ,d2J,Ce_ME,J_like,p] = calculateCost(signals,Ce,b,bp,Cb,T_c,config,dataLevel)
% Evaluate negative log-posterior, gradient, and Gauss-Newton Hessian.
% Combines parameter priors with data misfit from convolving the inferred
% impulse response with the input signal.
%
% Inputs:
%   signals   - Prepared signal struct from prepareSignals.
%   Ce        - Data-noise variance.
%   b         - Parameter vector in parameter space.
%   bp        - Prior mean vector in parameter space.
%   Cb        - Prior covariance matrix in parameter space.
%   T_c       - Reference time scale.
%   config    - Configuration struct.
%   dataLevel - Optional data level key ('coarse' or 'fine').
%
% Outputs:
%   J      - Total negative log-posterior.
%   dJ     - Gradient of J with respect to b.
%   d2J    - Gauss-Newton Hessian approximation.
%   Ce_ME  - Evidence-maximizing noise variance estimate.
%   J_like - Data likelihood contribution to J.
%   p      - Predicted output signal aligned with observed q.

    if nargin < 8
        dataLevel = 'coarse';
    end

    % Unpack signals
    u = signals.time_domain.(dataLevel).u;
    q = signals.time_domain.(dataLevel).q;
    valid = signals.time_domain.(dataLevel).valid;
    q = q(valid); % Truncate q to the valid region for discrete convolution
    b_th = signals.time_domain.(dataLevel).t_h ./ T_c;
    dt = signals.time_domain.(dataLevel).dt;
    
    % Extract the LFL settings
    LFL = config.model.LFL;
    LFL_sigma = config.model.LFL_sigma;
    isLFL = ~isempty(LFL);

    % Compute sizes
    Nd = length(q);
    Nb = length(b);
    P  = Nb/3;

    % Initialize cost and gradient
    J   = 0;
    dJ  = zeros(Nb,1);
    d2J = zeros(Nb,Nb);
    
    % -------------------------------------------------------------------------
    % Cost of the prior
    % -------------------------------------------------------------------------
    invCb = diag(1./diag(Cb)); % Cb is diagonal, cheap inverse    
    J = J + 0.5 * (b - bp)' * invCb * (b - bp);
    J = J + 0.5.*sum(log(2 .* pi .* diag(Cb))); % Contribution of normalizing constant
    dJ = dJ + invCb * (b - bp);

    % Prior over the LFL
    if isLFL
        nIdx = (1:3:Nb);
        n    = b(nIdx);
        one  = ones(P,1);

        err = (one'*n - LFL);
        s2  = LFL_sigma^2;

        J = J + 0.5 * (err)^2 / s2 + 0.5 * log(2 * pi * s2);

        dJ(nIdx) = dJ(nIdx) + (err / s2) * one;

        d2J(nIdx,nIdx) = d2J(nIdx,nIdx) + (1/s2) * (one*one');
    end

    
    % -------------------------------------------------------------------------
    % Cost of the data
    % -------------------------------------------------------------------------
    % Compute the predicted HRR
    [h, dhdb] = calculateImpulseResponse(b, b_th, T_c);

    p = conv(u, h, 'valid') .* dt;
    dpdb = nan(Nd,Nb);
    for nb = 1:Nb
        dpdb(:,nb) = conv(u, dhdb(:,nb),'valid') .* dt;
    end

    % Compute the residual discrepancy
    r = p - q;
    
    % Data likelihood: 0.5 * r' * inv(Ce) * r
    J_like = 0.5 * (r' * r) ./ Ce; 
    % Include the Gaussian normalization constant
    J_like = J_like + 0.5 .* Nd .* log(2 .* pi .* Ce); 

    J = J + J_like;

    % Compute gradient of cost w.r.t hrr prediction
    pJpp = r ./ Ce; 

    % Compute gradient of cost w.r.t parameters due to residual term
    pJpb_r = dpdb' * pJpp;

    % Add gradient of cost w.r.t parameters
    dJ = dJ + pJpb_r;

    % First order approximation of Hessian
    d2J = d2J + invCb + (dpdb' * dpdb) ./ Ce;

    % -------------------------------------------------------------------------
    % Compute the measurement noise that maximizes evidence, as per MacKay
    % 'Comparison of approximate methods for handling hyperparameters',
    % Neural Comp. 1999
    % -------------------------------------------------------------------------
    % Compute the effective number of parameters
    [cholR,flag] = chol(d2J);
    if flag > 0
        % Hessian is not positive definite, return effective number of
        % parameters equals actual number of parameters
        gamma = Nb;
    else
        X = cholR \ (cholR' \ invCb); % inv(d2J) * inv(Cb)
        gamma = Nb - trace(X); % Effective number of parameters
    end
    Ce_ME = (r' * r) ./ (Nd - gamma);
end