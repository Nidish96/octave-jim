clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')
addpath('../ROUTINES/export_fig/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

%% Parameters
M = eye(2);
K = [2 -1;-1 2];
C = 0.01*K;
b = 0.5;

MDL = MDOFGEN(M, K, C, eye(2));
MDL = MDL.SETNLFUN(1+3, [1 0], @(t, u, ud) deal(b*u.^3, 3*b*u.^2, ud*0));

[V, Wr] = eig(K, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

%%
rng(1);
Rk = randn(2,2);
resfun = @(l) deal(eigs(K+Rk-l*M, 1, 'SM'), -V(:,1)'*M*V(:,1));

l0 = 2;

opt = struct('Nmax', 100, 'Display', 1);
ls = NSOLVE(resfun, l0, opt)

%%
rfun = @(l) eigs(K-l*M, 1, 'SM');

rng(1)
lr = randn(1);

[r0, dr0] = resfun(lr);

hms = logspace(log10(2.649), log10(2.649015), 100);
drnum = arrayfun(@(hm) (rfun(lr+hm)-rfun(lr-hm))/2/hm, hms);

clf()
loglog(hms, abs(dr0-drnum), '.-')
grid on