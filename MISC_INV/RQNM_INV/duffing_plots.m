clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)
%% Load Data
load('./DATA/Duffing_EPMC.mat', 'UwxC');
load('./DATA/Duffing_RQNM.mat', 'Qs', 'Zts', 'Lams', 'Phi', 'GM', 'UlC');
[prev.Qp, si] = sort(UlC(end, UlC(end,:)>=0));
prev.Wp = sqrt(UlC(end-1, UlC(end,:)>=0));  prev.Wp = prev.Wp(si);
[prev.Qm, si] = sort(UlC(end, UlC(end,:)<=0));
prev.Wm = sqrt(UlC(end-1, UlC(end,:)<=0));  prev.Wm = prev.Wm(si);
load('./DATA/Duffing_FRESP.mat', 'UCs', 'Fl', 'Fas', 'Wst', 'Wen', 'h', 'Nhc');

%% Just Backbones
figure(1);
clf();
set(gcf, 'Color', 'white')

aa1 = plot(10.^UwxC(end,:), UwxC(end-2,:), 'k-', 'LineWidth', 2); hold on
aa2 = plot(Qs, sqrt(Lams), '-', 'LineWidth', 1);
aa3 = plot(prev.Qp, prev.Wp, '--', 'LineWidth', 2);
aa4 = plot(abs(prev.Qm), prev.Wm, '--', 'LineWidth', 2);

xlim([1e-2 1e2])
set(gca, 'Xscale', 'log')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (rad/s)')

% IS 1: [1.2,2.4] x [8200,8380]
% IS 2: [28,44] x [1.25,1.34]

plot([1.2, 1.2, 2.4, 2.4, 1.2], [8200, 8380, 8380, 8200, 8200], 'k-')
plot([28, 28, 44, 44, 28], [1.25, 1.34, 1.34, 1.25, 1.25]*1e4, 'k-')

set(gca, 'YTick', [0.8 1.1 1.4]*1e4)
legend([aa1, aa2, aa3, aa4], 'EPMC', 'RQNM', '$\sqrt{\lambda^+}$', '$\sqrt{\lambda^-}$', ...
    'Location', 'best', 'interpreter', 'latex')

is1 = axes('position', [0.22, 0.2, 0.25, 0.25]);
plot(10.^UwxC(end,:), UwxC(end-2,:), 'k-', 'LineWidth', 2); hold on
plot(Qs, sqrt(Lams), '-');
plot(prev.Qp, prev.Wp, '--', 'LineWidth', 2)
plot(abs(prev.Qm), prev.Wm, '--', 'LineWidth', 2)
xlim([1.2, 2.4]);
ylim([8200,8380]);
set(gca, 'XTick', [1.4 2.2])
set(gca, 'YTick', [8250 8350])

is2 = axes('position', [0.5, 0.6, 0.25, 0.25]);
plot(10.^UwxC(end,:), UwxC(end-2,:), 'k-', 'LineWidth', 2); hold on
plot(Qs, sqrt(Lams), '-');
plot(prev.Qp, prev.Wp, '--', 'LineWidth', 2)
plot(abs(prev.Qm), prev.Wm, '--', 'LineWidth', 2)
xlim([28, 44]);
ylim([1.25e4,1.34e4]);
set(gca, 'XTick', [30 40])
set(gca, 'YTick', [1.26e4 1.32e4])

figure(2);
clf()
plot(10.^UwxC(end,:), UwxC(end-1,:)./(2*UwxC(end-2,:))); hold on
plot(Qs, Zts, '-');

xlim([1e-2 1e2])
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'linear')
xlabel('Modal Amplitude')
ylabel('Damping Factor')

%% Against FRF
figure(3)
clf()

aa = gobjects(length(Fas), 1);
for fi=1:length(Fas)
    aa(fi) = plot(UCs{fi}(end,:), ...
        sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UCs{fi}(1:end-1, :)).^2)), ...
        '-', 'LineWidth', 2);
    hold on
    legend(aa(fi), sprintf('F = %.2f MN', Fas(fi)*1e-6));
end
aa1 = plot(UwxC(end-2,:), (10.^UwxC(end,:)).*sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UwxC(1:end-3, :)).^2)), ...
    'k-', 'LineWidth', 1.5); hold on
aa2 = plot(sqrt(Lams), sqrt(0.5)*abs(GM.NLTs.L*Phi).*Qs', '--', 'Color', [1 1 1]*0.4,...
    'LineWidth', 1.5);

legend(aa1, 'EPMC')
legend(aa2, 'RQNM')

legend([aa(1:end); aa1; aa2], 'Location', 'best')

set(gca, 'yscale', 'linear')
xlim([Wst Wen])
xlabel('Frequency (rad/s)')
ylabel('RMS Amplitude')

% export_fig('./SDOFFIGS/1_DUFFFRESP.eps', '-depsc')