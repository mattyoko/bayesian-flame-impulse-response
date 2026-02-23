function [hSI,FTFSI] = TFDSI(u,q,t,N,w)
% Estimate impulse response and transfer function using SI toolbox methods.
% Uses impulseest with regularization, then extracts impulse and FTF data.
%
% Inputs:
%   u - Input signal.
%   q - Output signal.
%   t - Time vector for impulse-response evaluation.
%   N - Model order / FIR horizon used by impulseest.
%   w - Angular frequencies for bode evaluation.
%
% Outputs:
%   hSI   - Struct with impulse estimate (val, time, std).
%   FTFSI - Struct with gain, phase, frequency, and uncertainty metrics.

dt = mean(diff(t));
id = iddata(q,u,dt);

impulseEstOptions = impulseestOptions;
impulseEstOptions.RegularizationKernel = 'CS'; % ['TC','none','CS','SE','SS','HF','DI','DC']
impulseEstOptions.PW = 'auto';
impulseEstOptions.Advanced.SearchMethod = 'fmincon';
impulseEstOptions.Advanced.AROrder = 5;
sys = impulseest(id,N,0,impulseEstOptions);

[hSI.val,hSI.time,~,hSI.std] = impulse(sys,t);
[FTFSI.gain,FTFSI.phase,FTFSI.w,FTFSI.std_gain,FTFSI.std_phase] = bode(sys,w);
FTFSI.gain = squeeze(FTFSI.gain);
FTFSI.phase = (squeeze(FTFSI.phase))*pi/180 - max(squeeze(FTFSI.phase))*pi/180;
FTFSI.std_gain = squeeze(FTFSI.std_gain);
FTFSI.std_phase = squeeze(FTFSI.std_phase)*pi/180;
end