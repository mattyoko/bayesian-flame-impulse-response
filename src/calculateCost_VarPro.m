function [J,dJ,d2J,Ce_ME,J_like,p,b_out] = calculateCost_VarPro(signals,Ce,x,bp,Cb,T_c,config,dataLevel)
% Evaluate the variable-projection objective for nonlinear parameters.
% Optimizes only x = [gamma1; beta1; ...], while linear amplitudes n are
% solved analytically at each call using a MAP linear solve.
%
% Inputs:
%   signals   - Prepared signal struct from prepareSignals.
%   Ce        - Data-noise variance.
%   x         - Nonlinear parameter vector [2N x 1].
%   bp        - Prior mean vector for full parameters.
%   Cb        - Prior covariance matrix for full parameters.
%   T_c       - Reference time scale.
%   config    - Configuration struct.
%   dataLevel - Optional data level key ('coarse' or 'fine').
%
% Outputs:
%   J      - Total negative log-posterior at the projected solution.
%   dJ     - Gradient with respect to x.
%   d2J    - Gauss-Newton Hessian approximation with respect to x.
%   Ce_ME  - Evidence-maximizing noise variance estimate.
%   J_like - Data likelihood contribution to J.
%   p      - Predicted output signal aligned with observed q.
%   b_out  - Full parameter vector with projected linear amplitudes.

    if nargin < 8
        dataLevel = 'coarse';
    end

    % Unpack signals
    u = signals.time_domain.(dataLevel).u;
    q = signals.time_domain.(dataLevel).q;
    valid = signals.time_domain.(dataLevel).valid;
    q = q(valid);
    b_th = signals.time_domain.(dataLevel).t_h ./ T_c;
    dt = signals.time_domain.(dataLevel).dt;

    % Extract the LFL settings
    LFL = config.model.LFL;
    LFL_sigma = config.model.LFL_sigma;
    isLFL = ~isempty(LFL);

    % Compute sizes and indices
    Nb          = length(bp);
    Nx          = Nb/3*2;
    Nd          = length(q);
    P           = Nb/3;
    ind_n_b     = 1:3:Nb;
    ind_beta_b  = 2:3:Nb;
    ind_gamma_b = 3:3:Nb;
    ind_beta_x  = 1:2:Nx;
    ind_gamma_x = 2:2:Nx;

    % Split the prior into linear, n, and nonlinear, x, parts
    np = bp(ind_n_b); Cn_diag = diag(Cb(ind_n_b,ind_n_b));
    xp = zeros(Nx,1); Cx_diag = zeros(Nx,1);
    xp(ind_beta_x) = bp(ind_beta_b);
    xp(ind_gamma_x) = bp(ind_gamma_b);
    Cx_diag(ind_beta_x) = diag(Cb(ind_beta_b,ind_beta_b));
    Cx_diag(ind_gamma_x) = diag(Cb(ind_gamma_b,ind_gamma_b));

    invCb = diag(1 ./ diag(Cb));
    invCx = diag(1 ./ Cx_diag);
    invCn = diag(1 ./ Cn_diag);

    % ---------------------------------------------------------------------
    % Solve for n_MAP given x
    % ---------------------------------------------------------------------
    b_tmp = zeros(Nb,1);
    b_tmp(ind_gamma_b) = x(ind_gamma_x); 
    b_tmp(ind_beta_b)  = x(ind_beta_x);
    [~, ~, G] = calculateImpulseResponse(b_tmp, b_th, T_c);

    % Build design matrix A such that p = A*n
    A = zeros(Nd,P);
    for ii = 1:P
        A(:,ii) = conv(u, G(:,ii), 'valid') .* dt;
    end

    % Solve linear system for n_MAP given 
    Hn = (A' * A) ./ Ce + invCn;
    gn = (A' * q) ./ Ce + invCn * np;
    if isLFL
        % Include a prior on sum(n)
        Hn = Hn + (1/LFL_sigma^2);
        gn = gn + (LFL/LFL_sigma^2);       
    end
    n = Hn \ gn;

    % Add linear terms to the full parameter vector, b_out 
    b_out = b_tmp;            
    b_out(ind_n_b) = n;  

    % ---------------------------------------------------------------------
    % Calculate the cost function at n_MAP, and its sensitivities w.r.t x
    % ---------------------------------------------------------------------
    % Initialize the cost, gradient and Hessian
    J   = 0;
    dJ  = zeros(Nx,1);
    d2J = zeros(Nx,Nx);

    % ------------------
    % Cost of the prior:
    % ------------------
    J = J + 0.5 * (b_out - bp)' * invCb * (b_out - bp);
    J = J + 0.5 * sum(log(2*pi*diag(Cb))); % Contribution of normalizing constant
    dJ = dJ + invCx * (x - xp);

    % Prior over the LFL
    if isLFL
        err = sum(n) - LFL;
        s2 = LFL_sigma^2;

        J = J + 0.5 * err^2 / s2 + 0.5 * log(2*pi*s2);

        % Gradient and Hessian contributions are zero
    end

    % ------------------
    % Cost of the data:
    % ------------------

    % Compute impulse response and sensitivities wrt b_out
    [h, dhdb] = calculateImpulseResponse(b_out, b_th, T_c);
    % Compute HRR prediction and sensitivities wrt x
    p = conv(u, h, 'valid') .* dt;
    dpdx = zeros(Nd,2*P);
    for ii = 1:P
        col_gamma = 3*ii - 1;
        col_beta  = 3*ii;
        dpdx(:,2*ii-1) = conv(u, dhdb(:,col_gamma), 'valid') .* dt;
        dpdx(:,2*ii  ) = conv(u, dhdb(:,col_beta), 'valid') .* dt;
    end
    % Compute the residual discrepancy
    r = p - q;

    % Data likelihood: 0.5 * r' * inv(Ce) * r
    J_like = 0.5 * (r' * r) ./ Ce + 0.5 * Nd * log(2*pi*Ce);
    J = J + J_like;

    % Compute gradient of cost w.r.t hrr prediction
    pJpp = r ./ Ce; 

    % Compute gradient of cost w.r.t parameters due to residual term
    pJpx = dpdx' * pJpp;
    dJ   = dJ + pJpx;

    % First order approximation of the Hessian
    d2J = d2J + invCx + (dpdx' * dpdx) ./ Ce;

    % -------------------------------------------------------------------------
    % Compute the measurement noise that maximizes evidence, as per MacKay
    % 'Comparison of approximate methods for handling hyperparameters',
    % Neural Comp. 1999
    % -------------------------------------------------------------------------
    % Compute the effective number of nonlinear parameters
    [cholR,flag] = chol(d2J);
    if flag > 0
        % Hessian is not positive definite, return effective number of
        % parameters equals actual number of parameters
        gamma_x = Nx;
    else
        X = cholR \ (cholR' \ invCx); % inv(d2J) * inv(Cb)
        gamma_x = Nx - trace(X); % Effective number of parameters
    end

    % Compute the effective number of linear parameters
    [cholR,flag] = chol(Hn);
    if flag > 0
        gamma_n = P;
    else
        X = cholR \ (cholR' \ invCn);
        gamma_n = P - trace(X);
    end

    % Total effective parameters
    gamma = gamma_x + gamma_n;

    Ce_ME = (r' * r) ./ (Nd - gamma);
end