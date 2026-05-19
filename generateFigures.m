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

%% !! Run this first !!
% Set paths, system properties and load the full dataset
ptd  = 'data/BRS_EderSilva23/';

% System properties:
L_ref = 50e-3; % Flame length
V_ref = 11.3; % Bulk flow velocity
T_c  = L_ref/V_ref; % Reference time scale for non-dimensionalization

% Load data
load(sprintf('%sdata_raw_incomp.mat',ptd));
exp_FTF = load(sprintf('%sFTF_exp_30kW_Front_lambda1.3.mat',ptd));

u_full = data_raw_incomp.u;
q_full = data_raw_incomp.y;
dt = data_raw_incomp.Ts;
t_full = (0:length(u_full)-1)*dt;
t_full_ms = t_full*1e3;

u_full = (u_full - mean(u_full))/mean(u_full);
q_full = (q_full - mean(q_full))/mean(q_full);

%% Figure 2
% Check git version is correct
checkGit();

% Generate default prior, as described in §4.1
config = loadDefaultConfig();
[mu_b,Cb] = generatePrior(1,1,config);

% Prepare plot
f2 = figure(2);
set(f2,'Units','normalized','Position',[0.1 0.1 0.6 0.3])
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

%% Figure 4
% Check git version is correct
checkGit();

% Prepare the figure
f5 = figure(4);
set(f5,'Units','normalized','Position',[0.1 0.1 0.5 0.4])

% Plot the normalized fluctuations
hu = plot(t_full_ms,u_full,'color',getColours(1),'lineWidth',1);
hold on
hq = plot(t_full_ms,q_full,'color',getColours(2),'lineWidth',1);
xlim([100, 130]);
grid on
box on
xlabel('time [ms]')
ylabel('normalized fluctuation')
legend([hu,hq],{'$u''/ \bar{u}$','$q''/ \bar{q}$'},'interpreter','latex','box','off');

%% Figures 5 - 7    
% Check git version is correct
checkGit();

% General settings
config = loadDefaultConfig();
config.preproc.DSmode = 'rate';
config.preproc.DSvalue = 3e-4;

% -------------------------------------------------------------------------
% Run Bayesian inference for baseline case (LFL free)
% -------------------------------------------------------------------------

[h, posterior] = inferImpulseResponse(u_full, q_full, t_full, T_c, modelOrders=1:5, config=config);

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
[h, posterior, modelRanking] = inferImpulseResponse(u_full, q_full, t_full, T_c, modelOrders=1:5, config=config);

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

signals = prepareSignals(u_full,q_full,1/dt,0,config);
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
f5 = plotModelComparisonH(modelRanking.logML,modelRanking.logOF,modelRanking.logBFL,...
    'modelLabels',modelLabels,'figUnits','Normalized','figWidth',0.6,'figHeight',0.5,'figNum',5);

f6=figure(6);
set(f6,'Units','normalized','Position',[0.1 0.1 0.5 0.4])
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

f7 = figure(7);
set(f7,'Units','normalized','Position',[0.1 0.1 0.5 0.6])
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

%% Figures 8 - 9
% Check git version is correct
checkGit();

% General settings
config = loadDefaultConfig();
config.preproc.DSmode = 'rate';
config.preproc.DSvalue = 3e-4;
config.model.T_h = 0.015;
config.model.LFL = 1;

fact = [1 0.2 0.15 0.05];
cols = getColours([1 4],length(fact));
colExp = getColours(2);

N_h = fix(config.model.T_h / dt);

% Initialize error table
factLabels = compose('%d%%',fact(2:end)*100);
N = length(fact)-1;
errMat = nan(N,4);

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

    % Compute FTF
    f = linspace(0,500);
    omega = 2*pi*f;
    FTF = calculateFTF(a,omega,Ca);

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
    
    f8=figure(8);
    set(f8,'Units','normalized','Position',[0.1 0.1 0.5 0.6])
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
    
    f9 = figure(9);
    set(f9,'Units','normalized','Position',[0.1 0.1 0.7 0.6])
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

    if ll > 1
        errMat(ll-1,1) = norm(FTFSI.gain - gainRef)/norm(gainRef);
        errMat(ll-1,2) = norm(FTF.gain - gainRef)/norm(gainRef);
        errMat(ll-1,3) = norm(FTFSI.phase - phaseRef)/norm(phaseRef);
        errMat(ll-1,4) = norm(FTF.phase - phaseRef)/norm(phaseRef);
    end
end

err = array2table(errMat,"RowNames",factLabels,"VariableNames",{'GainSI','GainBI','PhaseSI','PhaseBI'});
disp(err);

figure(8)
nexttile(1)
legend(hS1,compose('%d%%',fact*100),'box','off');

figure(9)
nexttile(2)
legend([hE, hB],['exp' compose('%d%%',fact*100)],'box','off');

%% Figures 10 - 11
% Check git version is correct
checkGit();

% General settings
config = loadDefaultConfig();
config.preproc.DSmode = 'rate';
config.preproc.DSvalue = 3e-4;
config.model.T_h = 0.015;
config.model.LFL = 1;
config.inference.runMCMC = true;
config.inference.MCMCiter = 500000;

% Set the data reduction factors
fact = [1 0.2 0.05];
cols = getColours(3,length(fact));
colExp = getColours(2);

% Set some constants
N_h = fix(config.model.T_h / dt);
N = 3;
D = 3*N;               
nbins = 31;

% Prepare the figure
f10 = figure(10);
set(f10,'Units','normalized','Position',[0.1 0.1 0.33 0.4])
t1 = tiledlayout(N, 3, 'TileSpacing','none', 'Padding','compact');

% Common z-grid and unit Gaussian (same for all tiles)
zEdges = linspace(-4, 4, nbins+1);
zAx    = linspace(-4, 4, 250);
phi0   = normpdf(zAx, 0, 1);
yMax   = 1.10*max(phi0);

% Parameter names
names = {'$n$','$\gamma$','$\beta$'};

for ll = 1:length(fact)
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
        figure(f10);
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

% Plot the full posterior for the 15% case
fact = 0.15;
M = N_h + fix((length(u_full) - N_h)*fact);

u = u_full(1:M);
q = q_full(1:M);
t = t_full(1:M);

config.inference.plotMCMC = true;

inferImpulseResponse(u, q, t, T_c, modelOrders=N, config=config);



%% Utilities
function checkGit()
    % Check that the current git version matches the required tag for paper-figure reproducibility.
    requiredTag = "paper-released";

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
