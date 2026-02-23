function FTF = calculateFTF(a,omega,Ca)
% Evaluate the flame transfer function from pulse-model parameters.
% Computes the MAP transfer function and, when covariance is provided,
% estimates 95% credible intervals by Monte Carlo sampling.
%
% Inputs:
%   a     - Physical parameter vector [3N x 1].
%   omega - Angular frequency vector in rad/s.
%   Ca    - Optional parameter covariance matrix [3N x 3N].
%
% Outputs:
%   FTF - Struct containing gain/phase and optional uncertainty bounds.

omega = omega(:);

P = length(a)/3;

FTF_cmplx = zeros(size(omega));

if nargin == 3
    N = 5000;
    aa = mvnrnd(a,Ca,N);
    aa = permute(aa,[2,3,1]);
    % Stack the MAP value as the first page of aa
    aa = cat(3,a,aa);
else
    aa = a;
end

for p = 1:P
    col_n   = 3*p - 2;
    col_tau = 3*p - 1;
    col_sig = 3*p;

    n   = aa(col_n,:,:);
    tau = aa(col_tau,:,:);
    sig = aa(col_sig,:,:);

    g = exp(-1i.*omega.*tau - 0.5.*(omega).^2.*sig.^2);

    FTF_cmplx = FTF_cmplx + n .* g;
end

gain = abs(FTF_cmplx);
phase = angle(FTF_cmplx);

% Extract the MAP point FTF
FTF.gain = gain(:,:,1);
phase0 = phase(:,:,1);
FTF.phase = unwrap(phase0);

if nargin == 3
    % Extract the ±95% gain credible interval centered around the MAP
    dgain = gain - FTF.gain;
    FTF.gain95lo = -quantile(dgain,0.025,3);
    FTF.gain95hi = +quantile(dgain,0.975,3);

    % Extract the ±95% phase credible interval centered around the MAP
    dphase = angle(exp(1i*(phase - phase0)));
    dphase_lo = quantile(dphase,0.025,3);
    dphase_hi = quantile(dphase,0.975,3);
    FTF.phase95lo = -dphase_lo;
    FTF.phase95hi = +dphase_hi;

    % Unwrap the phase
    FTF.phase95lo = FTF.phase - unwrap(FTF.phase - FTF.phase95lo);
    FTF.phase95hi = unwrap(FTF.phase + FTF.phase95hi) - FTF.phase;
end

