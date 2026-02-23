function [a,map,grad, Ca] = mapParameterToPhysicalSpace(b, t_scale, Cb)
% Map parameter-space variables into physical-space pulse parameters.
% Uses cumulative log-gap mapping for delays and exponential mapping for
% widths to enforce positivity.
%
% Inputs:
%   b       - Parameter vector in parameter space [3N x 1] 
%   t_scale - Reference time scale.
%   Cb      - Optional covariance in parameter space [3N x 3N].
%
% Outputs:
%   a    - Physical-space parameter vector [3N x 1] 
%   map  - Struct with fields a_n, a_tau, and a_sig.
%   grad - Struct with Jacobian blocks for each mapped component.
%   Ca   - Physical-space covariance (returned only when Cb is supplied).

% Compute sizes
Nb = size(b,1);
Ndel = Nb/3; % Number of delays
S = size(b,2);

% Unpack the parameters
b_n = b(1:3:end,:); b_gamma = b(2:3:end,:); b_beta = b(3:3:end,:);

% Map to physical space
map.a_n   = b_n;
a_gaps    = exp(b_gamma).*t_scale; 
map.a_sig = exp(b_beta).*t_scale;

% Convert a_gaps to a_tau
map.a_tau = cumsum(a_gaps);

% Pack back into interleaved parameter vector
a = cat(3,map.a_n,map.a_tau,map.a_sig);
a = reshape(permute(a,[3,1,2]),Nb,S);

% ----- Compute derivatives -----
% d(n)/d(n)
grad.a_n.pb_n = eye(Ndel);
% d(tau)/d(gamma)
grad.a_tau.pb_gamma = tril(exp(repmat(b_gamma',Ndel,1)).*t_scale);
% d(sigma)/d(beta)
grad.a_sig.pb_beta = diag(exp(b_beta).*t_scale);

% ----- Compute covariance in physical space if Cb provided -----
if nargin == 3
    % Build full Jacobian J_full = da / db
    J_full = zeros(Nb, Nb);

    idx_n   = (1:3:Nb).';
    idx_gam = (2:3:Nb).';
    idx_bet = (3:3:Nb).';

    % a_n rows
    J_full(idx_n, idx_n) = grad.a_n.pb_n;
    % a_tau rows 
    J_full(idx_gam, idx_gam) = grad.a_tau.pb_gamma;
    % a_sig rows
    J_full(idx_bet, idx_bet) = grad.a_sig.pb_beta;

    Ca = J_full * Cb * J_full';
end