function [h, dhdb, G, hVar] = calculateImpulseResponse(b, b_th, T_c, Cb)
% Compute impulse response values and parameter sensitivities for a pulse model.
% Maps parameters to physical space, evaluates the Gaussian pulse basis,
% and returns both h(t) and dh/db.
%
% Inputs:
%   b     - [3N x 1] parameter vector packed as
%           [n1; gamma1; beta1; n2; gamma2; beta2; ...].
%   b_th  - [T x 1] time vector in non-dimensional parameter space.
%   T_c   - Reference time scale for non-dimensionalization.
%   Cb    - Optional parameter covariance matrix [3N x 3N].
%
% Outputs:
%   h     - [T x 1] impulse response.
%   dhdb  - [T x 3N] Jacobian of h with respect to b.
%   G     - [T x N] Gaussian basis matrix in physical time.
%   hVar  - [T x 1] variance of h when Cb is provided.

% Unpack the parameters
b = b(:);
b_th = b_th(:);

% Map to physical space
[~,map,map_grad] = mapParameterToPhysicalSpace(b,T_c);
a_n   = map.a_n;
a_tau = map.a_tau;
a_sig = map.a_sig;
a_th  = b_th.*T_c;

% Determine some lengths
P = length(a_n);
T_h = length(b_th);

% Broadcast t - tau and sig
dt = a_th - a_tau'; 
S = a_sig';         

% Compute Gaussian G 
G = (1 ./ sqrt(2 * pi * S.^2)) .* exp(- (dt.^2) ./ (2 * S.^2));  % [T x N]

% Compute impulse response
h = G * a_n;  

% Compute gradients in physical space
dh_dan   = G;                                           
dh_datau = G .* (dt ./ (S.^2)) .* a_n';                 
dh_dasig = G .* (((dt.^2 - S.^2) ./ (S.^3))) .* a_n';   

% Convert gradients to parameter space using chain rule
dh_dbn      = dh_dan * map_grad.a_n.pb_n;              
dh_dbgamma  = dh_datau * map_grad.a_tau.pb_gamma; 
dh_dbbeta   = (dh_dasig * map_grad.a_sig.pb_beta);     

% Stack columns in parameter order [bn1, bgamma1, bbeta1, bn2, ...]
dhdb = zeros(T_h, 3*P);
for jj = 1:P
    cols = (3*jj-2):(3*jj);
    dhdb(:, cols) = [dh_dbn(:,jj), dh_dbgamma(:,jj), dh_dbbeta(:,jj)];
end

if nargin == 4
    % Cheaply compute variance of h, which is the diagonal of dhdb * Cb * dhdb'
    hVar = sum((dhdb * Cb) .* dhdb, 2); 
end