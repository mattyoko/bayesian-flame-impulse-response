function signals = prepareSignals(u, q, fs, T_h, config)
% Prepare coarse and fine signal representations for inference.
% Estimates a bandwidth from the input signal energy, then down-samples
% with anti-alias filtering where appropriate.
%
% Inputs:
%   u      - Input fluctuation signal (column vector).
%   q      - Output fluctuation signal (column vector).
%   fs     - Original sampling frequency in Hz.
%   T_h    - Desired impulse-response duration in seconds.
%   config - Configuration struct (uses config.preproc.DSlimit).
%
% Outputs:
%   signals - Struct with:
%             * signals.time_domain.coarse: possibly down-sampled data.
%             * signals.time_domain.fine: original-resolution data.

% Determine the integer downsampling factor
switch config.preproc.DSmode
    case 'factor'
        ds_factor = fix(config.preproc.DSvalue);
    case 'frequency'
        ds_factor = fix(fs ./ config.preproc.DSvalue);
    case 'rate'
        ds_factor = fix(fs .* config.preproc.DSvalue);
end

if ds_factor <= 1
    % No downsampling needed
    signals.time_domain.coarse = packageOutputs(u, q, fs, ds_factor, T_h);
    signals.time_domain.fine = signals.time_domain.coarse;
    return;
end

% Downsample signals with anti-aliasing
u_ds = safeResample(u, ds_factor);
q_ds = safeResample(q, ds_factor);
fs_ds = fs / ds_factor;

% Package outputs
signals.time_domain.coarse = packageOutputs(u_ds, q_ds, fs_ds, ds_factor, T_h);
signals.time_domain.fine = packageOutputs(u, q, fs, 1, T_h);



% !!! Legacy code from v1.0 - compute the PSD and downsample to just above
% Nyquist. May want to bring back at some point.

% % Compute single-sided spectrum of real-valued input signal, padding with
% % zeros to the nearest power of 2
% n_h = ceil(T_h * fs) + 1;  % desired impulse response length in samples
% n = length(u);
% n_padded = n + n_h - 1;
% Nfft = 2^nextpow2(n_padded);  % length of FFT after zero-padding
% 
% [U, f] = rfft(u, fs, Nfft);
% 
% % Estimate the cutoff frequency where 99.9% of the energy is retained
% psd = abs(U).^2;
% psd_cum = cumsum(psd);
% psd_cum = psd_cum / psd_cum(end);  % Normalize to [0,1]
% f_cut_idx = find(psd_cum >= 0.999, 1, 'first');
% if isempty(f_cut_idx)
%     f_cut = f(end);
% else
%     f_cut = f(f_cut_idx);
% end
% 
% % Add 10% margin to cutoff frequency
% f_cut = f_cut * 1.1;
% 
% % Compute integer downsampling factor (largest integer to satisfy Nyquist)
% ds_factor = max(1, floor(fs / (2 * f_cut))); 
% if ~isempty(config.preproc.DSlimit)
%     ds_factor = min(ds_factor,config.preproc.DSlimit);
% end
% 
% if ds_factor == 1
%     % No downsampling needed
%     signals.time_domain.coarse = packageOutputs(u, q, fs, ds_factor, T_h);
%     signals.time_domain.fine = signals.time_domain.coarse;
%     return;
% end

end

function signals = packageOutputs(u, q, fs, ds_factor, T_h)
    % Package time-domain signals and derived vectors into one struct.
    n   = length(u);
    n_h = ceil(T_h * fs) + 1;
    signals.u = u;
    signals.q = q;
    signals.fs = fs;
    signals.dt = 1 / fs;
    signals.t = (0:n-1)' * signals.dt;
    signals.t_h = (0:n_h-1)' * signals.dt;
    signals.n = n;
    signals.valid = n_h:n;
    signals.ds_factor = ds_factor;
end

function x_rs = safeResample(x,rs_factor)
    % Resample with symmetric edge padding to reduce boundary artifacts.

    P = 100;
    Q = P*rs_factor;
    x = x(:);
    N = length(x);
    % Compute new signal length
    l = floor(N*P/Q)+1;
    idx = l:ceil(N*P/Q)+l-1;
    % Embed the signal to avoid border effects of the resample function
    xx = [repmat(x(1,:),N,1); x ; repmat(x(end,:),N,1)];
    x_rs = resample(xx,P,Q);
    x_rs = x_rs(idx);
end

function [S, f, Nfft] = rfft(s, fs, Nfft)
    % Compute a single-sided FFT for a real-valued time signal.
    % Returns the retained non-negative frequency bins only.
    %
    % Inputs:
    %   s    - Real-valued time signal.
    %   fs   - Sampling frequency in Hz.
    %   Nfft - Optional FFT length used for zero-padding.
    %
    % Outputs:
    %   S    - Complex single-sided FFT bins.
    %   f    - Frequency vector corresponding to S.
    %   Nfft - FFT length used.
    
    if nargin < 3
        Nfft = length(s);
    end
    n = length(s);
    s_p = [s(:) ; zeros(Nfft - n, 1)];  % Zero-pad signal to length Nfft
    X = fft(s_p);          

    % find the indices for single-sided spectrum
    if rem(Nfft,2) == 0
        idx = (0:(Nfft/2))';            % includes Nyquist bin
    else
        idx = (0:floor(Nfft/2))';       % up to floor(Nfft/2)
    end

    S = X(idx+1);                       % single-sided bins only
    f  = fs * idx / Nfft;               % frequency vector
end