function [bMAP,CbMAP,hMAP,Ce,logML,logOF,logBFL] = estimatePosterior(signals,bp,Cb,T_c,config)
% Estimate MAP posterior quantities using multi-start nonlinear
% optimization. Reconstructs the full MAP parameter vector, and computes
% Laplace-based covariance and model scores. Optimizer is based on
% Nielsen's implementation of the Levenberg-Marquardt trust-region method
% [1], with modifications for variable projection and data noise
% estimation.
% 
% [1] Hans Bruun Nielsen, '﻿Damping parameter in Marquardt's method', IMM,
% 1999
% 
% Multiple starts can optionally be run in parallel (recommended if the
% parallel computing toolbox is available).
%
% Inputs:
%   signals - Prepared signal struct from prepareSignals.
%   bp      - Prior mean vector in parameter space.
%   Cb      - Prior covariance matrix in parameter space.
%   T_c     - Reference time scale.
%   config  - Configuration struct controlling optimizer and inference.
%
% Outputs:
%   bMAP   - MAP parameter vector in parameter space.
%   CbMAP  - Approximate posterior covariance at MAP.
%   hMAP   - Impulse-response struct (time, value, variance).
%   Ce     - Selected/estimated data-noise variance.
%   logML  - Log model evidence approximation.
%   logOF  - Objective-function-based score.
%   logBFL - Bayes free-energy lower bound approximation.

    % Unpack config
    rsrtScaling = config.optimizer.restartScaling;
    itrScaling  = config.optimizer.iterationScaling;
    tolOpt      = config.optimizer.optimalityTolerance;
    tolStep     = config.optimizer.stepSizeTolerance;
    tolFn       = config.optimizer.functionChangeTolerance;
    useParallel = config.optimizer.useParallel;  
    needLogML   = strcmp(config.optimizer.solnSelection,'logML');
    dataVariance= config.inference.dataVariance;
    inferNoise  = isempty(dataVariance);

    % Get the size of the full parameter vector,b, then calculate the
    % number of nonlinear parameters, x
    Nb_full = numel(bp);
    Nx = Nb_full/3*2;
    P  = numel(bp) / 3;
    
    % Settings
    maxIter = itrScaling * Nx; % Maximum number of iterations
    rstrts  = rsrtScaling * Nx;% Number of random restarts
    lamInc0 = 2;               % Base lambda increment for rejected step
    lamDec  = 3;               % Base lambda decrement for accepted step
    p       = 3;               % Lambda scaling exponent  

    % Set initial uncertainty
    if inferNoise
        Ce0 = 1e-4;
    else
        Ce0 = dataVariance;
    end

    % Optimize only nonlinear parameters x = [gamma1;beta1;...]

    % Build xp,Cx for restart seeds
    [xp, Cx_diag] = extractNonlinPrior(bp,Cb,P);
    Cx = diag(Cx_diag);
    I  = eye(Nx);

    % Generate restart seeds over x
    x_samples = generateRestartSeeds(xp,Cx,rstrts);

    % -------------------------
    % Preallocate outputs
    % -------------------------
    J_hist      = inf(1,rstrts);
    Ce_hist     = nan(1,rstrts);
    J_like_hist = nan(1,rstrts);
    logML_hist  = -inf(1,rstrts);
    b_hist      = nan(Nb_full, rstrts);
    BFGS_correction = nan(Nx,Nx,rstrts);

    % -------------------------
    % Parallel (or serial) restart loop
    % -------------------------
    if useParallel
        parfor sp = 1:rstrts
            [J_hist(sp),Ce_hist(sp),J_like_hist(sp),logML_hist(sp),b_hist(:,sp),BFGS_correction(:,:,sp)] = ...
                runOneRestart(sp, x_samples(:,sp), Ce0, bp, Cb, T_c, config, ...
                              maxIter, tolOpt, tolStep, tolFn, lamInc0, lamDec, p, I, needLogML, inferNoise, signals);
        end
    else
        for sp = 1:rstrts
            [J_hist(sp),Ce_hist(sp),J_like_hist(sp),logML_hist(sp),b_hist(:,sp),BFGS_correction(:,:,sp)] = ...
                runOneRestart(sp, x_samples(:,sp), Ce0, bp, Cb, T_c, config, ...
                              maxIter, tolOpt, tolStep, tolFn, lamInc0, lamDec, p, I, needLogML, inferNoise, signals);
        end
    end

    % Select best
    switch config.optimizer.solnSelection
        case 'logProb'
            [~,indMin] = min(J_hist);
        case 'logML'
            [~,indMin] = max(logML_hist);
    end

    bMAP = b_hist(:,indMin);
    Ce   = Ce_hist(indMin);
    J_MAP = J_hist(indMin);
    J_like = J_like_hist(indMin);
    BFGS_correction = BFGS_correction(:,:,indMin);

    % Select the Hessian to be used for Laplace and model ranking
    switch config.optimizer.finalHessianMethod
        case 'GN'
            % Calculate the full Gauss-Newton hessian with one call to
            % calculateCost without variable projection
            [~,~,Hfull] = calculateCost(signals,Ce,bMAP,bp,Cb,T_c,config);
        case 'BFGS'
            % Calculate the full Gauss-Newton hessian with one call to
            % calculateCost without variable projection
            [~,~,Hfull] = calculateCost(signals,Ce,bMAP,bp,Cb,T_c,config);
            % Find the entries of Hfull that correspond to the nonlinear
            % parameters, x
            ind_x_in_b = reshape([3*(1:P)-1; 3*(1:P)], [], 1); % [gamma1; beta1; gamma2; beta2; ...]
            % Apply the BFGS curvature correction to the nonlinear entries
            % (linear entries are exact from Gauss-Newton)
            Hfull(ind_x_in_b, ind_x_in_b) = Hfull(ind_x_in_b, ind_x_in_b) + BFGS_correction;
        otherwise
            Hfull = FDHessian(signals, Ce, bMAP, bp, Cb, T_c, config);
    end

    % Robust inverse of Hessian to get posterior covariance
    Hfull = (Hfull + Hfull')/2;
    [R,pc] = chol(Hfull);
    if pc > 0
        warning('Hessian is not positive definite. Adding jitter for numerical stability')
        jitter = 1e-6 * eye(length(bMAP));
        for ii = 1:10
            [R,pc] = chol(Hfull + jitter);
            if pc == 0, break; end
            jitter = 10 * jitter;
        end
        if pc > 0
            error('Hessian is not positive definite even after adding jitter. Cannot compute posterior covariance')
        end
    end
    CbMAP = R \ (R' \ eye(length(bMAP)));

    % Calculate the impulse response and its uncertainty
    b_th = signals.time_domain.coarse.t_h ./ T_c;
    b_th_eval = linspace(0,b_th(end),500);
    [hMAP.val, ~, ~, hMAP.var]  = calculateImpulseResponse(bMAP, b_th_eval, T_c, CbMAP);
    hMAP.time = b_th_eval;

    logML  = -J_MAP + 0.5*length(bMAP)*log(2*pi) - sum(log(diag(R)));
    logBFL = -J_like;
    logOF  = logML - logBFL;
end

function [J,Ce,J_like,logML,b_full,BFGS_correction,exitMessage] = runOneRestart(~, x0, Ce0, bp, Cb, T_c, config, ...
                                                                               maxIter, tolOpt, tolStep, tolFn, lamInc0, lamDec, p, I, needLogML, inferNoise, signals)
    Ce = Ce0;
    x  = x0;

    % Initial cost
    [f, g, H, ~, J_like, ~, b_full] = calculateCost_VarPro(signals,Ce,x,bp,Cb,T_c,config);

    % BFGS curvature approximation
    A = 1e-10*eye(size(H));

    % Initial damping
    lambda = max(diag(H));
    lamInc = lamInc0;

    db = ones(size(x));
    df = 1;

    exitMessage = '';

    for k = 1:maxIter
        if k > 5
            if norm(g) < tolOpt
                exitMessage = sprintf('Converged at iteration %d because the optimality condition was met', k);
                break;
            end
            if (norm(db)/max(norm(x),eps)) < tolStep
                exitMessage = sprintf('Converged at iteration %d because the step size reduced below the tolerance', k);
                break;
            end
            if norm(df) < tolFn
                exitMessage = sprintf('Converged at iteration %d because the function change dropped below the tolerance', k);
                break;
            end
        end

        H_reg = H + lambda * I;
        condH = cond(H_reg);
        if (condH > 1e20) || isnan(condH)
            exitMessage = 'Restart failed due to ill-conditioned Hessian';
            break;
        end

        db = -H_reg \ g;
        predictedImprovement = 0.5 * db' * (lambda * db - g);

        [f_new, g_new, H_new, Ce_new, J_like_new, ~, b_full_new] = calculateCost_VarPro(signals, Ce, x + db, bp, Cb, T_c, config);

        actualImprovement = f - f_new;

        if predictedImprovement <= 0 
            rho = -inf;
        else
            rho = actualImprovement / predictedImprovement;
        end

        if rho > 0 % Accept step
            
            x_old = x; g_old = g;

            x      = x + db;
            b_full = b_full_new;

            df = f - f_new;
            f  = f_new; g = g_new; H = H_new;
            H  = (H + H') / 2;
            J_like = J_like_new;

            % Scale trust region 
            scale  = max(1./lamDec, 1 - (lamInc0 - 1)*(2*rho - 1)^p);
            lambda = max(lambda * scale,1e-8);
            lamInc = lamInc0;

            % BFGS approximation to the second order Hessian component
            if norm(g) < 1
                s = x - x_old;
                y = g - g_old;
                yGN = H * s;
                u = y - yGN;

                aa = 1 ./ (u' * s);
                bb = 1 ./ (s' * A * s);
                v  = A * s;
                if s'*u > 1e-12
                    A = A + aa * (u * u') - bb * (v * v');
                end
            end

            % Step towards MML data noise
            if inferNoise && (k > 5)
                Ce = 0.6*Ce + 0.4*Ce_new;
                [f, g, H, ~, J_like, ~, b_full] = calculateCost_VarPro(signals, Ce, x, bp, Cb, T_c, config);
                H = (H + H')/2;
            end
        else % Reject step
            % Scale trust region
            lambda = lambda * lamInc;
            lamInc = 2 * lamInc;
        end
    end

    if isempty(exitMessage)
        if k == maxIter
            exitMessage = 'Maximum iterations reached';
        else
            exitMessage = 'Terminated';
        end
    end

    % Final outputs
    J     = f;
    BFGS_correction = A;

    % Optional per-restart logML (needed if selecting by logML)
    if needLogML
        J_full = calculateCost(signals, Ce, b_full, bp, Cb, T_c, config);
        H_full = FDHessian(signals, Ce, b_full, bp, Cb, T_c, config);
        H_full = (H_full + H_full')/2;
        [R,pc] = chol(H_full);
        if pc == 0
            logML = -J_full + 0.5*length(b_full)*log(2*pi) - sum(log(diag(R)));
        else
            logML = -inf;
        end
    else
        logML = -inf;
    end
end

function [xp, Cx_diag] = extractNonlinPrior(bp,Cb,P)
    % Extract the nonlinear parameters, x, from the full prior vector, b.

    Cb_diag = diag(Cb);
    xp = zeros(2*P,1);
    Cx_diag = zeros(2*P,1);

    for j = 1:P
        xp(2*j-1) = bp(3*j-1);
        xp(2*j  ) = bp(3*j  );
        Cx_diag(2*j-1) = Cb_diag(3*j-1);
        Cx_diag(2*j  ) = Cb_diag(3*j  );
    end
end

function H = FDHessian(signals, Ce, b, bp, Cb, T_c, config, dataLevel)

    if nargin < 8
        dataLevel = 'coarse';
    end

    % Build log-posterior at fixed Ce
    logPost = @(x) -calculateCost(signals, Ce, x, bp, Cb, T_c, config, dataLevel);

    % Numerical Hessian at MAP
    Nb = length(b);
    H = zeros(Nb);
    eps = 1e-5; 
    for ii = 1:Nb
        ei = zeros(Nb,1); ei(ii) = 1;
        for jj = ii:Nb
            ej = zeros(Nb,1); ej(jj) = 1;
            if ii == jj
                lp_p = logPost(b + eps*ei);
                lp_m = logPost(b - eps*ei);
                lp_0 = logPost(b);
                H(ii,ii) = - (lp_p - 2*lp_0 + lp_m) / (eps^2);
            else
                lp_pp = logPost(b + eps*ei + eps*ej);
                lp_pm = logPost(b + eps*ei - eps*ej);
                lp_mp = logPost(b - eps*ei + eps*ej);
                lp_mm = logPost(b - eps*ei - eps*ej);
                H(ii,jj) = - (lp_pp - lp_pm - lp_mp + lp_mm) / (4*eps^2);
                H(jj,ii) = H(ii,jj);
            end
        end
    end
end