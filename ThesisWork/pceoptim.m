clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/QUADRATURE/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

pdists = [makedist('exp') repmat(makedist('normal'), 1, 6)];
polfuns = [{@(ns, xs) PLAGU(ns, xs)}, repmat({@(ns, xs) PHERM(ns, xs)}, 1, 6)];

mdi = 1; 
zeta0 = 1.3841e-4;

inds = [1 2 3 4 6 7];
pref = sprintf('nlbb_%s', sprintf('%d', inds));
Lx = eye(7);  Lx(:, 5) = [];
%%
EXPdat = load(sprintf('./MATFILES/Mode%d_High.mat', mdi), 'AMP_avg', 'FRE_avg', 'DAM_avg');
EXPdat.AMP_avg = EXPdat.AMP_avg./(2*pi*EXPdat.FRE_avg).^2;
EXPdat.FRE_avg = 2*pi*EXPdat.FRE_avg;
ni = find(isfinite(EXPdat.AMP_avg));
EXPdat.AMP_avg = EXPdat.AMP_avg(ni); 
EXPdat.FRE_avg = EXPdat.FRE_avg(ni);
EXPdat.DAM_avg = EXPdat.DAM_avg(ni);

fname = sprintf('./ALLPCE/%s_cofs.mat', pref);
load(fname, 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'IJs', 'Integs')
Zcofs(:,1) = Zcofs(:,1)+zeta0;
%%  COST FUNCTION
[dw] = PCECOSTFUN([zeros(1, 6); 0.5*ones(1, 6); ones(1, 6)], Lx, {Rcofs, Wcofs, Zcofs}, EXPdat, IJs, polfuns);

%% gamultiobj
Npop = 500;
xs = zeros(Npop, 7);
for i=1:7
    xs(:, i) = random(pdists(i), Npop, 1);
end
xs = xs(:, inds);

opt = optimoptions('gamultiobj', 'PopulationSize', Npop, 'UseVectorized', true, 'UseParallel', true, 'InitialPopulationMatrix', xs, 'PlotFcn', {'gaplotpareto', 'gaplotparetodistance', 'gaplotrankhist', 'gaplotspread'});
[x, fval, ~, op] = gamultiobj(@(x) PCECOSTFUN(x, Lx, {Rcofs, Wcofs, Zcofs}, EXPdat, IJs, polfuns), ...
    6, [], [], [], [], [0 repmat(-nan, 1, 5)], [], [], opt);

%% Plot Results
figure(3); 
clf(); 
for i=1:size(x,1)
    subplot(2,1,1); 
    semilogx(EVALPCE(x(i,:)*Lx', Rcofs, IJs, polfuns), EVALPCE(x(i,:)*Lx', Wcofs/2/pi, IJs, polfuns), '-');
    hold on
    
    subplot(2,1,2); 
    semilogx(EVALPCE(x(i,:)*Lx', Rcofs, IJs, polfuns), EVALPCE(x(i,:)*Lx', Zcofs, IJs, polfuns), '-');
    hold on
end
subplot(2,1,1); 
semilogx(EXPdat.AMP_avg, EXPdat.FRE_avg/2/pi, 'k-', 'LineWidth', 2);  hold on
subplot(2,1,2); 
semilogx(EXPdat.AMP_avg, EXPdat.DAM_avg, 'k-', 'LineWidth', 2); hold on

%% 
figure(4)
clf()
plotmatrix(x)

%% Conduct Original Simulations
xrun = x(1,:);
[~, mi] = min(x(:,1));
xrun = x(mi, :);
pref = 'GADES';
RQNM_EXPRSURF_PCEFUN(xrun*Lx', 0, ones(1, 7), pref, 1, [-7.5, -3], 1, 'proper');

%% Plot
model = 'BRB_Thesis';
load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'R');
load(sprintf('./ALLPCE/%s_%d_m%d.mat', pref, 0, mdi), 'Qs', 'Phi', 'Lams', 'Zts', 'Wstat', 'Ustat');

%%
function [wzres] = PCECOSTFUN(x, Lx, RWZcofs, EXPdat, IJs, polfuns)
    
    if size(x, 2)~=(size(IJs,2)-1)
        x = x';
    end
    x = x*Lx';
%     Rs = EVALPCE(x, RWZcofs{1}, IJs, polfuns);
%     Ws = EVALPCE(x, RWZcofs{2}, IJs, polfuns);
%     Zs = EVALPCE(x, RWZcofs{3}, IJs, polfuns);
    
    Psi = 1;  % PCE Bases
    for j=1:size(IJs, 2)
        psi = polfuns{j}(IJs(:,j), x(:, j));
        Psi = Psi.*psi;
    end
    Rs = RWZcofs{1}*Psi';
    Ws = RWZcofs{2}*Psi';
    Zs = RWZcofs{3}*Psi';
    
    wres = zeros(size(x,1),1);
    zres = zeros(size(x,1),1);
    for i=1:size(x,1)
        Wp = interp1(Rs(:,i), Ws(:,i), EXPdat.AMP_avg);
        Zp = interp1(Rs(:,i), Zs(:,i), EXPdat.AMP_avg);
        
        si = find(isfinite(Wp) & isfinite(Zp));
    
        wres(i) = rms((Wp(si)-EXPdat.FRE_avg(si))./EXPdat.FRE_avg(si));
        zres(i) = rms(Zp(si)-EXPdat.DAM_avg(si));
    end
    
    wzres = [wres, zres];
end