clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

mds = [1 3 5];
model = 'BRB_Thesis';
zoomlims = [175 181; 591 596; 1195 1202];
Pamps = [1e-7 2e-5 1e-4; 1e-8 5e-6 1e-4; 1e-9 5e-6 1e-4];

Wstat_th = zeros(3, 1);
Wstat_exp = zeros(3, 1);

zeta0 = 1.3841e-4;
Zt_exp = zeros(3, 1);
%% Mesh Structure
Nds = dlmread(sprintf('../MODELS/%s/Nodes.dat', model));
Quad = dlmread(sprintf('../MODELS/%s/Elements.dat', model));

Nq = 2;
MESH = MESH2D(Nds, 3, [], Quad, Nq);

%% Plots
for mdi=1:3
exp(1) = load(sprintf('./MATFILES/Mode%d_Low.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');
exp(2) = load(sprintf('./MATFILES/Mode%d_Med.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');
exp(3) = load(sprintf('./MATFILES/Mode%d_High.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');

load(sprintf('./ALLPCE/MEANMODELS/mmbb_0_m%d.mat', mdi), 'Qs', 'Lams', ...
    'Zts', 'Phi', 'Wstat', 'Tstat', 'Dfluxes');
% load(sprintf('./ALLPCE/FLATMEANMODELS/mmbb_0_m%d.mat', mdi), 'Qs', 'Lams', ...
%     'Zts', 'Phi', 'Wstat');
load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'R');

%% Static Traction
if mdi==1
    for i=1:3
        figure((mdi-1)*10+i)
        clf()
        set(gcf, 'Color', 'white')
        MESH.SHOWFIELD2D(Tstat(i:3:end))
        axis equal
        xx=colorbar('south');
        xlabel(xx, sprintf('Traction %d (N m$^{-2}$)', i), 'interpreter', 'latex')
        pos = get(xx, 'Position');
        pos(2) = 0.24;
        set(xx, 'Position', pos)
        axis off
        colormap(jet)
%         export_fig(sprintf('./FIGS/MM_STATTRAC_%d.png', i), '-dpng');
        export_fig(sprintf('./FIGS/MM_STATTRAC_%d.eps', i), '-depsc');
    end
end

%%
Rs = (abs(R(3,:)*Phi)'.*Qs);
[pis, ni] = find((Rs(1:end-1)-Pamps(mdi,:)).*(Rs(2:end)-Pamps(mdi,:))<=0);
if length(pis)==2
    pis = [pis; length(Rs)];
end
%% Backbones
figure((mdi-1)*10+4)
clf()
set(gcf, 'Color', 'white')
% set(gcf, 'Position', [2800 550 1200 480])
% subplot(1,2, 1)
for i=1:3
    semilogx(exp(i).AMP_avg./(2*pi*exp(i).FRE_avg).^2, exp(i).FRE_avg, 'k.'); hold on
end
semilogx(Rs, sqrt(Lams)/2/pi, 'b-', 'LineWidth', 2);
semilogx(Rs(pis), sqrt(Lams(pis))/2/pi, 'ro', 'MarkerFaceColor', 'w', 'LineWidth', 2);
xlim([min(Rs) min(1e-4,max(Rs))])

xlabel('Response Amplitude (m)')
ylabel('Natural Frequency (Hz)')
grid on

aax = axes('Position', [0.20 0.20 0.36 0.3]);
for i=1:3
    semilogx(exp(i).AMP_avg./(2*pi*exp(i).FRE_avg).^2, exp(i).FRE_avg, 'k.'); hold on
end
semilogx(Rs, sqrt(Lams)/2/pi, 'b-', 'LineWidth', 2);
semilogx(Rs(pis), sqrt(Lams(pis))/2/pi, 'ro', 'MarkerFaceColor', 'w', 'LineWidth', 2);
xlim([min(Rs) min(1e-4,max(Rs))])
ylim(zoomlims(mdi,:))
% xlim([min(Rs) max(Rs)])
export_fig(sprintf('./FIGS/MM_BBW_mode%d.eps', mdi), '-depsc')

figure((mdi-1)*10+5)
clf()
set(gcf, 'Color', 'white')
% subplot(1,2, 2)
for i=1:3
    aae = semilogx(exp(i).AMP_avg./(2*pi*exp(i).FRE_avg).^2, exp(i).DAM_avg*100, 'k.'); hold on
end
aam = semilogx(Rs, (Zts+zeta0)*100, 'b-', 'LineWidth', 2);
semilogx(Rs(pis), (Zts(pis)+zeta0)*100, 'ro', 'MarkerFaceColor', 'w', 'LineWidth', 2);

legend([aae aam], 'Experimental Measurements', 'Model Predictions', ...
       'Location', 'best');
xlabel('Response Amplitude (m)')
ylabel('Damping Factor (\%)')
grid on
xlim([min(Rs) min(1e-4,max(Rs))])
% xlim([min(Rs) max(Rs)])
export_fig(sprintf('./FIGS/MM_BBZ_mode%d.eps', mdi), '-depsc')

% export_fig(sprintf('./MEANMODELRESP/BBFIG_M%d.png', mdi), '-dpng');

%% Dissipation Fluxes
% for pj=1:length(pis)
%     for i=1:3
%         figure((mdi-1)*10+5+pj)
%         clf()
%         set(gcf, 'Color', 'white')
%         MESH.SHOWFIELD2D(abs(MESH.Qm\Dfluxes(i:3:end, pis(pj)))*2*pi/sqrt(Lams(pis(pj))))
%         axis equal
%         xx=colorbar('south');
%         xlabel(xx, sprintf('Cyclic Dissipation Flux %d (J m$^{-2}$)', i), 'interpreter', 'latex')
%         pos = get(xx, 'Position');
%         pos(2) = 0.24;
%         set(xx, 'Position', pos)
%         axis off
%         colormap(jet)
%         set(gca, 'colorscale', 'linear')
%         export_fig(sprintf('./FIGS/MM_M%d_DFLUX_P%d_C%d.eps', mdi, pj, i), '-depsc');
%         
%         set(gca, 'colorscale', 'log')
%         export_fig(sprintf('./FIGS/MM_M%d_DFLUXL_P%d_C%d.eps', mdi, pj, i), '-depsc');
%     end
% end

%% Linear (Low-Amplitude) Parameters
fav = robustfit([exp(1).AMP_avg; exp(2).AMP_avg; exp(3).AMP_avg], ...
                [exp(1).FRE_avg; exp(2).FRE_avg; exp(3).FRE_avg]);
zav = robustfit([exp(1).AMP_avg; exp(2).AMP_avg; exp(3).AMP_avg], ...
                [exp(1).DAM_avg; exp(2).DAM_avg; exp(3).DAM_avg]);
Wstat_exp(mdi) = fav(1);
Wstat_th(mdi) = Wstat(mds(mdi))/2/pi;

Zt_exp(mdi) = zav(1);

% subplot(1,2,1); plot(xlim, fav(1)*[1 1], 'r-');
% subplot(1,2,2); plot(xlim, zav(1)*100*[1 1], 'r-');
end

%% Table
table(Zt_exp, Wstat_exp, Wstat_th, (Wstat_exp-Wstat_th)./Wstat_exp*100)
