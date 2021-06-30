clc
% clear all

addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/CONTACTMODELS/')
addpath('../ROUTINES/QUASISTATIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

model = 'BRB_Thesis';
E = 1.9231e11;
nu = 0.3;

% model = 'BRB';
% E = 2e11;
% nu = 0.3;

top   = 'R05B_Before';
bot   = 'R05A_Before';

load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'K', 'M', 'R', 'L', 'Fv', 'TFM');

%% Load Mesh
Nds = dlmread(sprintf('../MODELS/%s/Nodes.dat', model));
Quad = dlmread(sprintf('../MODELS/%s/Elements.dat', model));

Nq = 2;
MESH = MESH2D(Nds, 3, [], Quad, Nq);

%% Prepare Contact Model Parameters
MESH = MESH.SETQUAD(1);
Aels = full(sum(MESH.Tm));  % Element Areas
Aint = sum(Aels);
Aels = kron(Aels(:), ones(Nq^2,1));
MESH = MESH.SETQUAD(Nq);

load(sprintf('./MATFILES/%s_R1_AspPDEs.mat', top), 'PS_sds');
R1top = load(sprintf('./MATFILES/%s_R1_AspPDEs.mat', top), 'BilinPlaneQPs', 'LLX0s_sd', 'CRAD', 'NASPS', 'PS_sds');
R2top = load(sprintf('./MATFILES/%s_R2_AspPDEs.mat', top), 'BilinPlaneQPs', 'LLX0s_sd', 'CRAD', 'NASPS');

R1bot = load(sprintf('./MATFILES/%s_R1_AspPDEs.mat', bot), 'BilinPlaneQPs', 'LLX0s_sd', 'CRAD', 'NASPS');
R2bot = load(sprintf('./MATFILES/%s_R2_AspPDEs.mat', bot), 'BilinPlaneQPs', 'LLX0s_sd', 'CRAD', 'NASPS');

%% Gap Function
gap1 = R1top.BilinPlaneQPs(:,1)-R1bot.BilinPlaneQPs(:,1);
gap2 = R2top.BilinPlaneQPs(:,1)-R2bot.BilinPlaneQPs(:,1);
gap1 = gap1-max(gap1);
gap2 = gap2-max(gap2);

%%
figure(1)
clf()
MESH.SHOWFIELD2D(gap1)
axis equal 
colormap(jet);
xx=colorbar('south');
xlabel(xx, 'Gap Function (m)', 'interpreter', 'latex')
axis off
pos = get(xx, 'Position');
pos(2) = pos(2)+0.1;
set(xx, 'Position', pos);
set(gcf, 'Color', 'white')
export_fig('./FIGS/GAPFUN_R1.png', '-dpng', '-r300')

figure(2)
clf()
MESH.SHOWFIELD2D(gap2)
axis equal 
colormap(jet);
xx=colorbar('south');
xlabel(xx, 'Gap Function (m)', 'interpreter', 'latex')
axis off
pos = get(xx, 'Position');
pos(2) = pos(2)+0.1;
set(xx, 'Position', pos);
set(gcf, 'Color', 'white')
export_fig('./FIGS/GAPFUN_R2.png', '-dpng', '-r300')

%% Lambda
lam1 = (R1top.NASPS(:,1)+R1bot.NASPS(:,1))./(R1top.NASPS(:,1)./R1top.LLX0s_sd(:,1)+R1bot.NASPS(:,1)./R1bot.LLX0s_sd(:,1));
lam2 = (R2top.NASPS(:,1)+R2bot.NASPS(:,1))./(R2top.NASPS(:,1)./R2top.LLX0s_sd(:,1)+R2bot.NASPS(:,1)./R2bot.LLX0s_sd(:,1));

%%
figure(3)
clf()
MESH.SHOWFIELD2D(lam1)
axis equal 
colormap(jet);
xx=colorbar('south');
xlabel(xx, 'Asperity Height Exponent ($m^{-1}$)', 'interpreter', 'latex')
axis off
pos = get(xx, 'Position');
pos(2) = pos(2)+0.1;
set(xx, 'Position', pos);
set(gcf, 'Color', 'white')
caxis([0 2e6])
export_fig('./FIGS/LAM_R1.png', '-dpng', '-r300')

figure(4)
clf()
MESH.SHOWFIELD2D(lam2)
axis equal 
colormap(jet);
xx=colorbar('south');
xlabel(xx, 'Asperity Height Exponent ($m^{-1}$)', 'interpreter', 'latex')
axis off
pos = get(xx, 'Position');
pos(2) = pos(2)+0.1;
set(xx, 'Position', pos);
caxis([0 2e6])
set(gcf, 'Color', 'white')
export_fig('./FIGS/LAM_R2.png', '-dpng', '-r300')

%% Curvature Radii
R1 = (R1top.CRAD(:,1).*R1top.NASPS(:,1)+R1bot.CRAD(:,1).*R1bot.NASPS(:,1))./(R1top.NASPS(:,1)+R1bot.NASPS(:,1));
R2 = (R2top.CRAD(:,1).*R2top.NASPS(:,1)+R2bot.CRAD(:,1).*R2bot.NASPS(:,1))./(R2top.NASPS(:,1)+R2bot.NASPS(:,1));

%%
figure(5)
clf()
MESH.SHOWFIELD2D(R1)
axis equal 
colormap(jet);
xx=colorbar('south');
xlabel(xx, 'Asperity Mean Radius (m)', 'interpreter', 'latex')
axis off
pos = get(xx, 'Position');
pos(2) = pos(2)+0.1;
set(xx, 'Position', pos);
set(gcf, 'Color', 'white')
set(gca, 'ColorScale', 'log')
caxis([1e-4 1e-1])
export_fig('./FIGS/CRAD_R1.png', '-dpng', '-r300')

figure(6)
clf()
MESH.SHOWFIELD2D(R2)
axis equal 
colormap(jet);
xx=colorbar('south');
xlabel(xx, 'Asperity Mean Radius (m)', 'interpreter', 'latex')
axis off
pos = get(xx, 'Position');
pos(2) = pos(2)+0.1;
set(gcf, 'Color', 'white')
set(xx, 'Position', pos);
set(gca, 'ColorScale', 'log')
caxis([1e-4 1e-1])
export_fig('./FIGS/CRAD_R2.png', '-dpng', '-r300')

%%
figure(7)
clf()
plot(gap1*1e3, gap2*1e3, '.'); hold on
grid on
% axis equal
xlim([-4e-1 0])
ylim([-4e-1 0])
set(gca, 'YTick', get(gca, 'XTick'))
xlabel('Gap Function 1 (mm)')
ylabel('Gap Function 2 (mm)')
set(gcf, 'Color', 'white')
olis = (isoutlier(gap1) | isoutlier(gap2));
hold on; plot(gap1(olis)*1e3, gap2(olis)*1e3, 'o')
[rho, pval] = corr(gap1(~olis),gap2(~olis));
title(sprintf("(Pearson $\\rho$, p-value): (%.4f, %.4f)", rho, pval))
export_fig('./FIGS/GAPCORR.png', '-dpng', '-r300')

figure(8)
clf()
plot(lam1*1e-6, lam2*1e-6, '.'); hold on
grid on
% axis equal
xlim([0 2.5])
ylim([0 2.5])
xlabel('Asperity Height Exponent 1 ($\mu m^{-1}$)')
ylabel('Asperity Height Exponent 2 ($\mu m^{-1}$)')
set(gcf, 'Color', 'white')
olis = (isoutlier(lam1) | isoutlier(lam2));
hold on; plot(lam1(olis)*1e-6, lam2(olis)*1e-6, 'o')
[rho, pval] = corr(lam1(~olis),lam2(~olis));
title(sprintf("(Pearson $\\rho$, p-value): (%.4f, %.4f)", rho, pval))
export_fig('./FIGS/LAMCORR.png', '-dpng', '-r300')

figure(9)
clf()
loglog(R1*1e3, R2*1e3, '.'); hold on
grid on
% axis equal
xlim([1e-1 2e2])
ylim([1e-1 2e2])
xlabel('Asperity Mean Radius 1 (mm)')
ylabel('Asperity Mean Radius 2 (mm)')
set(gcf, 'Color', 'white')
olis = (isoutlier(R1) | isoutlier(R2));
hold on; plot(R1(olis)*1e3, R2(olis)*1e3, 'o')
[rho, pval] = corr(R1(~olis),R2(~olis));
title(sprintf("(Pearson $\\rho$, p-value): (%.4f, %.4f)", rho, pval))
export_fig('./FIGS/RADCORR.png', '-dpng', '-r300')