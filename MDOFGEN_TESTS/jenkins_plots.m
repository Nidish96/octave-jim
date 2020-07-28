clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)
%% Load Data
load('./DATA/Jenkins_EPMC.mat', 'UwxC');
load('./DATA/Jenkins_RQNM.mat', 'Qs', 'Zts', 'Lams', 'Phi', 'GM');
load('./DATA/Jenkins_FRESP.mat', 'UCs', 'Fl', 'Fas', 'Wst', 'Wen', 'h', 'Nhc');

%% Just Backbones
figure(1);
clf();
aa1 = plot(10.^UwxC(end,:), UwxC(end-2,:)); hold on
aa2 = plot(Qs, sqrt(Lams), '-');

xlim([1e-2 1e2])
set(gca, 'Xscale', 'log')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (rad/s)')

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
ylim([0 3.5])
xlabel('Frequency (rad/s)')
ylabel('RMS Amplitude')
