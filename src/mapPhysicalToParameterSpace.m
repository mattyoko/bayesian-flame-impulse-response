function b = mapPhysicalToParameterSpace(a,t_scale)
% Map physical-space pulse parameters back to parameter space.
% Inverts cumulative delays into log-gap parameters and non-dimensionalizes
% pulse widths via logarithmic scaling.
%
% Inputs:
%   a       - Physical-space parameter vector [3N x 1]
%   t_scale - Reference time scale.
%
% Outputs:
%   b       - Parameter-space vector [3N x 1]

% Compute sizes
Na = size(a,1);
S = size(a,2);

% Unpack the parameters
a_n = a(1:3:end,:); a_tau = a(2:3:end,:); a_sig = a(3:3:end,:);

% Calculate delay gaps
gaps = diff([0 ; a_tau]);

b_n   = a_n;
b_alp = log(gaps./t_scale);
b_sig = log(a_sig./t_scale);

% Pack back into interleaved parameter vector
b = cat(3,b_n,b_alp,b_sig);
b = reshape(permute(b,[3,1,2]),Na,S);
end