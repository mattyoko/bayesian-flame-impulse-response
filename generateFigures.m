% Generate the figures from Yoko & Polifke, 'Bayesian inference of flame
% impulse responses', CnF [under review], 2026.
% 
% This demonstrates the Bayesian method on input-output data from LES on
% the BRS burner, as reported in Eder et al., '﻿Incompressible versus
% compressible large eddy simulation for the identification of premixed
% flame dynamics', International Journal of Spray and Combustion, 2023.
% 
% Before running, ensure that the src, utils and data folders are added to
% your path.

%% Set paths
ptd  = 'data/BRS_EderSilva23/';

%% Figure 1
% Check git version is correct
checkGit();

% Generate default prior, as described in §4.1
config = loadDefaultConfig();
config.prior.mu = [0.00, 0.00, -1.80];
config.prior.sig= [1.00, 0.50, +0.50];
[mu_b,Cb] = generatePrior(1,1,config);

% Prepare plot
f1 = figure(1);
set(f1,'Units','normalized','Position',[0.1 0.1 0.6 0.3])
tiledlayout(1,3,'padding','compact','TileSpacing','compact');

% Plot prior over n
mu_n = mu_b(1);
sig_n = sqrt(Cb(1,1));
n = linspace(-5,5,200);
pn = normpdf(n,mu_n,sig_n);

nexttile(1)
plot(n,pn,'color',getColours(1),'lineWidth',1)
xlabel('$n_i$',interpreter='latex')
ylabel('$p(n_i)$',interpreter='latex')
title('(a) Prior over $n$',interpreter='latex')

% Plot prior over alpha
gam = linspace(-10,log(6),300);
mu_gam = mu_b(2);
sig_gam = sqrt(Cb(2,2));
pgam = normpdf(gam,mu_gam,sig_gam);
alpha = exp(gam);

nexttile(2)
plot(alpha, pgam,'color',getColours(1),'lineWidth',1)
xlabel('${\alpha}_i$', interpreter='latex')
ylabel('$p({\alpha}_i)$', interpreter='latex')
title('(b) Prior over $\alpha$', interpreter='latex')

% Plot the prior over sigma/T_c
bet = linspace(-10,log(1),300);
mu_beta = mu_b(3);
sig_beta = sqrt(Cb(3,3));
pbet = normpdf(bet,mu_beta,sig_beta);
s = exp(bet);

nexttile(3)
plot(s, pbet,'color',getColours(1),'lineWidth',1)
xlabel('$\tilde{\sigma}_i$', interpreter='latex')
ylabel('$p(\tilde{\sigma}_i)$', interpreter='latex')
title('(c) Prior over $\tilde{\sigma}$', interpreter='latex')

%% Figure 3
% Check git version is correct
checkGit();

% Load and normalize the data
load(sprintf('%sdata_raw_incomp.mat',ptd));
u = data_raw_incomp.u;
q = data_raw_incomp.y;
dt = data_raw_incomp.Ts;
t = (0:length(u)-1)*dt;
t_ms = t*1e3; % Time in ms

u = (u - mean(u))/mean(u);
q = (q - mean(q))/mean(q);

% Prepare the figure
f3 = figure(3);
set(f3,'Units','normalized','Position',[0.1 0.1 0.5 0.4])

% Plot the normalized fluctuations
hu = plot(t_ms,u,'color',getColours(1),'lineWidth',1);
hold on
hq = plot(t_ms,q,'color',getColours(2),'lineWidth',1);
xlim([100, 130]);
grid on
box on
xlabel('time [ms]')
ylabel('normalized fluctuation')
legend([hu,hq],{'$u''/ \bar{u}$','$q''/ \bar{q}$'},'interpreter','latex','box','off');

%% Figures 4 - 6    
% Check git version is correct
checkGit();

% General settings
config = loadDefaultConfig();
config.prior.mu = [0.00, 0.00, -1.80];
config.prior.sig= [1.00, 0.50, +0.50];
config.preproc.DSlimit = 300;

% System properties:
L_ref = 50e-3; % Flame length
V_ref = 11.3;  % Bulk flow velocity
T_c   = L_ref/V_ref; % Reference time scale for non-dimensionalization

% Load and normalize data
load(sprintf('%sdata_raw_incomp.mat',ptd));
exp_FTF = load(sprintf('%sFTF_exp_30kW_Front_lambda1.3.mat',ptd));

u = data_raw_incomp.u;
q = data_raw_incomp.y;
dt = data_raw_incomp.Ts;
t = (0:length(u)-1)*dt;

u = (u - mean(u))/mean(u);
q = (q - mean(q))/mean(q);

% -------------------------------------------------------------------------
% Run Bayesian inference for baseline case (LFL free)
% -------------------------------------------------------------------------

[h, posterior] = inferImpulseResponse(u, q, t, T_c, modelOrders=1:5, config=config);

% Unpack posterior struct
a = posterior.a;
Ca = posterior.Ca;

% convert h.time to ms time and reduce its amplitude
hms = h; hms.time = hms.time*1e3;
hms.val = hms.val/1e3;
hms.var = hms.var/1e6;

% Calculate FTF
f = linspace(0,500);
omega = 2*pi*f;
FTF = calculateFTF(a,omega,Ca);

% Pack results for plotting later
baseline.hms = hms;
baseline.FTF = FTF;
baseline.a = a;
baseline.Ca = Ca;

% -------------------------------------------------------------------------
% Run Bayesian inference for LFL = 1
% -------------------------------------------------------------------------

% Update settings
config.model.LFL = 1; % Set low-frequency limit to 1

% Run Bayesian inference
[h, posterior, modelRanking] = inferImpulseResponse(u, q, t, T_c, modelOrders=1:5, config=config);

% Unpack posterior struct
a = posterior.a;
Ca = posterior.Ca;

% convert h.time to ms and reduce its amplitude
hms = h; hms.time = hms.time*1e3;
hms.val = hms.val/1e3;
hms.var = hms.var/1e6;

FTF = calculateFTF(a,omega,Ca);

% -------------------------------------------------------------------------
% Run SI using the settings from Eder et al. 2023
% -------------------------------------------------------------------------

signals = prepareSignals(u,q,1/dt,0,config);
uds = signals.time_domain.coarse.u;
qds = signals.time_domain.coarse.q;
tds = signals.time_domain.coarse.t;
N = 39;
[hSI,FTFSI] = TFDSI(uds,qds,tds,N,omega);

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
colBI = getColours(1);
colExp = getColours(2);
colBaseline = getColours(4);
colSI = getColours(3);

modelLabels = {'N=1','N=2','N=3','N=4','N=5'};
f4 = plotModelComparisonH(modelRanking.logML,modelRanking.logOF,modelRanking.logBFL,...
    'modelLabels',modelLabels,'figUnits','Normalized','figWidth',0.6,'figHeight',0.5,'figNum',4);

f5=figure(5);
set(f5,'Units','normalized','Position',[0.1 0.1 0.5 0.4])
padding = 0.1;
hold on
hS  = errorpatch(hSI.time*1e3,hSI.val/1e3,2*hSI.std/1e3,2*hSI.std/1e3,'color',colSI);
hbl = errorpatch(baseline.hms.time,baseline.hms.val,2*sqrt(baseline.hms.var),2*sqrt(baseline.hms.var),'color',colBaseline);
hB  = errorpatch(hms.time,hms.val,2*sqrt(hms.var),2*sqrt(hms.var),'color',colBI);
xlabel('time [ms]')
ylabel('$h\times 10^{-3}$','interpreter','latex')
legend([hS,hB,hbl],'SI','BI LFL = 1','BI LFL free','box','off');
xlim([0,15]);
grid on
box on
drawnow

f6 = figure(6);
set(f6,'Units','normalized','Position',[0.1 0.1 0.5 0.6])
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile(1); hold on
hE  = plot(exp_FTF.freq_exp,exp_FTF.gain_exp,'o','color',colExp);
hS  = errorpatch(FTFSI.w/2/pi,FTFSI.gain,2*FTFSI.std_gain,2*FTFSI.std_gain,'color',colSI);
hbl = errorpatch(f,baseline.FTF.gain,baseline.FTF.gain95lo,baseline.FTF.gain95hi,'color',colBaseline);
hB  = errorpatch(f,FTF.gain,FTF.gain95lo,FTF.gain95hi,'color',colBI);
ylabel('$|\mathcal{F}|$','interpreter','latex')
grid on
box on
legend([hE,hS,hB,hbl],'Exp','SI','BI LFL = 1','BI LFL free','box','off');
nexttile(2); hold on
plot(exp_FTF.freq_exp,exp_FTF.phase_exp,'o','color',colExp)
errorpatch(FTFSI.w/2/pi,unwrap(FTFSI.phase),2*FTFSI.std_phase,2*FTFSI.std_phase,'color',colSI);
errorpatch(f,baseline.FTF.phase,baseline.FTF.phase95lo,baseline.FTF.phase95hi,'color',colBaseline);
errorpatch(f,FTF.phase,FTF.phase95lo,FTF.phase95hi,'color',colBI);
xlabel('frequency [Hz]')
ylabel('$\angle \mathcal{F}$ [rad]','interpreter','latex')
grid on
box on
drawnow

%% Figures 7 - 8
% Check git version is correct
checkGit();

% General settings
config = loadDefaultConfig();
config.prior.mu = [0.00, 0.00, -1.80];
config.prior.sig= [1.00, 0.50, +0.50];
config.preproc.DSlimit = 300;
config.model.T_h = 0.015;
config.model.LFL = 1;

% System properties:
L_ref = 50e-3; % Flame length
V_ref = 11.3; % Bulk flow velocity
T_c  = L_ref/V_ref; % Reference time scale for non-dimensionalization

% Load data
load(sprintf('%sdata_raw_incomp.mat',ptd));
exp_FTF = load(sprintf('%sFTF_exp_30kW_Front_lambda1.3.mat',ptd));

u = data_raw_incomp.u;
q = data_raw_incomp.y;
dt = data_raw_incomp.Ts;
t = (0:length(u)-1)*dt;

u = (u - mean(u))/mean(u);
q = (q - mean(q))/mean(q);

u_full = u;
q_full = q;
t_full = t;

fact = [1 0.2 0.1 0.05];
cols = getColours([1 4],length(fact));
colExp = getColours(2);

N_h = fix(config.model.T_h / dt);

for ll = 1:length(fact)

    M = N_h + fix((length(u_full) - N_h)*fact(ll));
    
    u = u_full(1:M);
    q = q_full(1:M);
    t = t_full(1:M);
    
    % Run Bayesian inference
    [h, posterior] = inferImpulseResponse(u, q, t, T_c, modelOrders=3, config=config);
    
    % Unpack posterior
    a = posterior.a;
    Ca = posterior.Ca;

    % convert h.time to ms and reduce its amplitude
    hms = h; hms.time = hms.time*1e3;
    hms.val = hms.val/1e3;
    hms.var = hms.var/1e6;

    % Run SI
    signals = prepareSignals(u,q,1/dt,0,config);
    uds = signals.time_domain.coarse.u;
    qds = signals.time_domain.coarse.q;
    tds = signals.time_domain.coarse.t;
    dtds = mean(diff(tds));
    N = fix((config.model.T_h + 2*dtds)/dtds);
    f = linspace(0,500);
    omega = 2*pi*f;
    [hSI,FTFSI] = TFDSI(uds,qds,tds,N,omega);
    
    f7=figure(7);
    set(f7,'Units','normalized','Position',[0.1 0.1 0.5 0.6])
    if ll == 1
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
    end
    nexttile(1)
    hold on
    hS1(ll) = errorpatch(hSI.time*1e3,hSI.val/1e3,2*hSI.std/1e3,2*hSI.std/1e3,'color',cols(ll,:));
    ylabel('$h\times 10^{-3}$','interpreter','latex')
    xlim([0,15]);
    grid on
    box on

    nexttile(2)
    hold on
    errorpatch(hms.time,hms.val,2*sqrt(hms.var),2*sqrt(hms.var),'color',cols(ll,:));
    xlabel('time [ms]')
    ylabel('$h\times 10^{-3}$','interpreter','latex')
    xlim([0,15]);
    grid on
    box on
    
    f = linspace(0,500);
    omega = 2*pi*f;
    FTF = calculateFTF(a,omega,Ca);
    f8 = figure(8);
    set(f8,'Units','normalized','Position',[0.1 0.1 0.7 0.6])
    if ll == 1
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        nexttile(1); hold on
        plot(exp_FTF.freq_exp,exp_FTF.gain_exp,'o','color',colExp);
        nexttile(2); hold on
        hE = plot(exp_FTF.freq_exp,exp_FTF.gain_exp,'o','color',colExp);
        nexttile(3); hold on
        plot(exp_FTF.freq_exp,exp_FTF.phase_exp,'o','color',colExp);
        nexttile(4); hold on
        plot(exp_FTF.freq_exp,exp_FTF.phase_exp,'o','color',colExp);
    end
    nexttile(1)
    errorpatch(FTFSI.w/2/pi,FTFSI.gain,2*FTFSI.std_gain,2*FTFSI.std_gain,'color',cols(ll,:));
    ylabel('$|\mathcal{F}|$','interpreter','latex')
    grid on
    box on
    ylim([0 3])
    nexttile(2); hold on
    hB(ll) = errorpatch(f,FTF.gain,FTF.gain95lo,FTF.gain95hi,'color',cols(ll,:));
    % ylabel('gain')
    grid on
    box on
    ylim([0 3])
    ylabel('$|\mathcal{F}|$','interpreter','latex')
    nexttile(3)
    errorpatch(FTFSI.w/2/pi,unwrap(FTFSI.phase),2*FTFSI.std_phase,2*FTFSI.std_phase,'color',cols(ll,:));
    xlabel('frequency [Hz]')
    ylabel('$\angle{}\mathcal{F}$ [rad]','interpreter','latex')
    grid on
    box on
    nexttile(4); hold on
    errorpatch(f,FTF.phase,FTF.phase95lo,FTF.phase95hi,'color',cols(ll,:));
    xlabel('frequency [Hz]','interpreter','latex')
    ylabel('$\angle{}\mathcal{F}$ [rad]','interpreter','latex')
    grid on
    box on
    
end

figure(7)
nexttile(1)
legend(hS1,compose('%d%%',fact*100),'box','off');

figure(8)
nexttile(2)
legend([hE, hB],['exp' compose('%d%%',fact*100)],'box','off');

%% Figures 9 - 10
% Check git version is correct
checkGit();

% General settings
config = loadDefaultConfig();
config.prior.mu = [0.00, 0.00, -1.08];
config.prior.sig= [1.00, 0.50, +0.50];
config.preproc.DSlimit = 300;
config.model.T_h = 0.015;
config.model.LFL = 1;
config.inference.runMCMC = true;
config.inference.MCMCiter = 200000;

% System properties:
L_ref = 50e-3; % Flame length
V_ref = 11.3; % Bulk flow velocity
T_c  = L_ref/V_ref; % Reference time scale for non-dimensionalization

% Load data
load(sprintf('%sdata_raw_incomp.mat',ptd));
exp_FTF = load(sprintf('%sFTF_exp_30kW_Front_lambda1.3.mat',ptd));

u = data_raw_incomp.u;
q = data_raw_incomp.y;
dt = data_raw_incomp.Ts;
t = (0:length(u)-1)*dt;

u = (u - mean(u))/mean(u);
q = (q - mean(q))/mean(q);

u_full = u;
q_full = q;
t_full = t;

% Set the data reduction factors
fact = [0.2 0.1 0.05];
cols = getColours(3,length(fact));
colExp = getColours(2);

N_h = fix(config.model.T_h / dt);

N = 3;
D = 3*N;               
nbins = 30;

f9 = figure(1);
set(f9,'Units','normalized','Position',[0.1 0.1 0.33 0.4])
t1 = tiledlayout(N, 3, 'TileSpacing','none', 'Padding','compact');

% Common z-grid and unit Gaussian (same for all tiles)
zEdges = linspace(-4, 4, nbins+1);
zAx    = linspace(-4, 4, 250);
phi0   = normpdf(zAx, 0, 1);
yMax   = 1.10*max(phi0);

% Parameter names
names = {'$n$','$\gamma$','$\beta$'};

for ll = 1:length(fact)
    if ll == 2
        config.inference.plotMCMC = true;
    else
        config.inference.plotMCMC = false;
    end
    M = N_h + fix((length(u_full) - N_h)*fact(ll));

    u = u_full(1:M);
    q = q_full(1:M);
    t = t_full(1:M);

    % Run Bayesian inference
    N = 3;
    [h, posterior, ~, mcmc] = inferImpulseResponse(u, q, t, T_c, modelOrders=N, config=config);

    % Laplace posterior
    b  = posterior.b;
    Cb = posterior.Cb;

    % MCMC posterior
    samp = mcmc.postSamples';     % [Nsamp x D]

    % ---- Plot standardized marginals z = (x-mu)/sig ----
    for pp = 1:D
        figure(f9);
        ax = nexttile(pp);
        hold on

        mu  = b(pp);
        sig = sqrt(Cb(pp,pp));

        % Guard against degenerate/invalid sigmas
        if ~isfinite(sig) || sig <= 0
            sig = 1;
        end

        z = (samp(:,pp) - mu) ./ sig;

        histogram(z, zEdges, 'Normalization','pdf', ...
            'EdgeColor', 'none', ...
            'FaceColor', cols(ll,:),...
            'FaceAlpha', 0.5);

        histogram(z, zEdges, 'Normalization','pdf', ...
            'DisplayStyle','stairs', ...
            'EdgeColor', cols(ll,:), ...
            'EdgeAlpha', 1.0,...
            'LineWidth', 1.0);

        if ll == length(fact)
            plot(zAx, phi0, 'Color', getColours(1), 'LineWidth', 2.0);
        end

        xlim([-4 4]);
        ylim([0 yMax]);

        % Remove axes/ticks (keep only shapes)
        set(ax, 'XTick',[], 'YTick',[]);
        box on

        if any(pp == 1:3)
            title(names{pp},'interpreter','latex');
        end
        if any(pp == [1,4,7])
            i_val = find(pp == [1,4,7]);
            ylabel(sprintf('i = %d',i_val),'interpreter','latex')
        end
    end
end



%% Utilities
function checkGit()
    % Check that the current git version matches the required tag for paper-figure reproducibility.
    requiredTag = "paper-v1.0";

    % Fast: succeed only if HEAD is exactly at the tag (no history walk).
    cmd = sprintf('git -C "%s" describe --tags --exact-match', pwd);
    [st, out] = system(cmd);

    if st == 0 && string(strtrim(out)) == requiredTag
        return
    else
        error("generateFigures:NotAtFrozenTag", ...
            "For reproducability, generateFigures.m must be run from the tag '%s'. \n\nRun:\n  git checkout %s\n", ...
            requiredTag, requiredTag);
    end
end
