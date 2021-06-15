clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

mds = [1 3 5];
model = 'BRB_Thesis';
zoomlims = [175 181; 591 596; 1195 1202];

Wstat_th = zeros(3, 1);
Wstat_exp = zeros(3, 1);

Zt_exp = zeros(3, 1);
%% Plots
for mdi=1:3
exp(1) = load(sprintf('./MATFILES/Mode%d_Low.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');
exp(2) = load(sprintf('./MATFILES/Mode%d_Med.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');
exp(3) = load(sprintf('./MATFILES/Mode%d_High.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');

load(sprintf('./MEANMODELRESP/MM_MODE%d.mat', mdi), 'Qs', 'Lams', ...
     'Zts', 'Phi', 'Wstat');
load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'R');

Rs = (abs(R(3,:)*Phi)'.*Qs).*Lams/9.81;

figure((mdi-1)*1+1)
clf()
set(gcf, 'Color', 'white')
set(gcf, 'Position', [2800 550 1200 480])
subplot(1,2, 1)
for i=1:3
    semilogx(exp(i).AMP_avg, exp(i).FRE_avg, 'k.'); hold on
end
semilogx(Rs, sqrt(Lams)/2/pi, 'b-', 'LineWidth', 2);
xlim([min(Rs) max(Rs)])

xlabel('Response Amplitude (g)')
ylabel('Natural Frequency (Hz)')
grid on

aax = axes('Position', [0.18 0.20 0.18 0.3]);
for i=1:3
    semilogx(exp(i).AMP_avg, exp(i).FRE_avg, 'k.'); hold on
end
semilogx(Rs, sqrt(Lams)/2/pi, 'b-', 'LineWidth', 2);
ylim(zoomlims(mdi,:))
xlim([min(Rs) max(Rs)])

subplot(1,2, 2)
for i=1:3
    aae = semilogx(exp(i).AMP_avg, exp(i).DAM_avg*100, 'k.'); hold on
end
aam = semilogx(Rs, Zts*100, 'b-', 'LineWidth', 2);

legend([aae aam], 'Experimental Measurements', 'Model Predictions', ...
       'Location', 'best');
xlabel('Response Amplitude (g)')
ylabel('Damping Factor (\%)')
grid on
xlim([min(Rs) max(Rs)])

export_fig(sprintf('./MEANMODELRESP/BBFIG_M%d.png', mdi), '-dpng');

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
