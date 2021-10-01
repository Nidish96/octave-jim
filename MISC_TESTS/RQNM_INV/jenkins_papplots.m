clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

cfg = 1;

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',16)
%% Load Data
load(sprintf('./DATA/Jenkins_EPMC_cfg%d.mat', cfg), 'UwxC');
load(sprintf('./DATA/Jenkins_RQNM_cfg%d.mat', cfg), 'As', 'Cs', 'Omegas', 'Phi', 'GM', 'Qs', 'UlC');
Ls = cellfun(@(c) c(end-1,:), UlC, 'UniformOutput', false);
clear UlC
load(sprintf('./DATA/Jenkins_FRESP_cfg%d.mat', cfg), 'UCs', 'Fl', 'Fas', 'Wst', 'Wen', 'h', 'Nhc');

%% Just Backbones
for mi=1:2
    figure((mi-1)*2+1);
    clf();  
    aa1 = plot(10.^UwxC{mi}(end,:), UwxC{mi}(end-2,:), 'k-', 'LineWidth', 2); hold on
    aa2 = plot(Qs, Omegas(:,mi), 'b-', 'LineWidth', 2);
%     aa3 = plot(abs(As{mi}), sqrt(Ls{mi}))

    xlim([1e-2 1e2])
    set(gca, 'Xscale', 'log')
    xlabel('Modal Amplitude')
    ylabel('Natural Frequency (rad/s)')

	yl = ylim;
    set(gca, 'YTick', fix(yl(1):(yl(2)-yl(1))/4:yl(2)))
    
    if mi==1 && cfg==1
        plot([1e-1 1e0 1e0 1e-1 1e-1], [8255 8255 8325 8325 8255], 'k--')
    end
    
    if mi==1 && cfg==1
        legend([aa1, aa2], 'EPMC', 'RQNM', 'Location', 'best', 'fontsize', 16)
    end
    
    print(sprintf('./FIGS/NMA_WM%d_cfg%d.eps', mi, cfg), '-depsc')

    figure((mi-1)*2+2);
    clf()
    plot(10.^UwxC{mi}(end,:), UwxC{mi}(end-1,:)./(2*UwxC{mi}(end-2,:)), 'k-', 'LineWidth', 2); hold on
    plot(Qs, Cs(:, mi)./(2*Omegas(:, mi)), 'b-', 'LineWidth', 2);

    if mi==1 && cfg==1
        plot([1e-1 1e0 1e0 1e-1 1e-1], [12 12 16 16 12]*1e-3, 'k--')
    end
    
    xlim([1e-2 1e2])
    yl = ylim;
%     ylim([0 yl(2)])
    set(gca, 'YTick', yl(1):(yl(2)-yl(1))/4:yl(2))
    set(gca, 'Xscale', 'log')
    set(gca, 'Yscale', 'linear')
    xlabel('Modal Amplitude')
    ylabel('Damping Factor')
    print(sprintf('./FIGS/NMA_ZM%d_cfg%d.eps', mi, cfg), '-depsc')
end

%% Modal Coupling
if cfg==1
    figure(5)
    clf();
    aa1 = plot(10.^UwxC{1}(end,:), UwxC{1}(end-2,:), 'k-', 'LineWidth', 2); hold on
    aa2 = plot(10.^UwxC{2}(end,:), UwxC{2}(end-2,:)/3, 'r--', 'LineWidth', 2); hold on
    plot([1e-1 1e0 1e0 1e-1 1e-1], [8255 8255 8325 8325 8255], 'k--')
    
    xlabel('Modal Amplitude')
    ylabel('Natural Frequency (rad/s)')
    
    yl = ylim;
    set(gca, 'YTick', fix(yl(1):(yl(2)-yl(1))/4:yl(2)))
    set(gca, 'xscale', 'log')
    
    legend([aa1 aa2], 'Mode 1', 'Mode 2 (div. by 3)', 'fontsize', 16)
    print('./FIGS/NMA_WM12_cfg1.eps', '-depsc')
    
    figure(6)
    clf();
    aa1 = plot(10.^UwxC{1}(end,:), UwxC{1}(end-2,:), 'k-', 'LineWidth', 2); hold on
    aa2 = plot(10.^UwxC{2}(end,:), UwxC{2}(end-2,:)/3, 'r--', 'LineWidth', 2); hold on
    
    xlabel('Modal Amplitude')
    ylabel('Natural Frequency (rad/s)')
    
    xlim([1e-1 1e0])
    ylim([8255 8325])
    yl = ylim;
    set(gca, 'YTick', fix(yl(1):(yl(2)-yl(1))/4:yl(2)))
    set(gca, 'xscale', 'log')
    print('./FIGS/NMA_ZM12_cfg1.eps', '-depsc')
end
%% Against FRF
figure(7)
clf()

aa = gobjects(length(Fas), 1);
for fi=1:length(Fas)
    aa(fi) = plot(UCs{fi}(end,:), ...
        sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UCs{fi}(1:end-1, :)).^2)), ...
        '-', 'LineWidth', 2);
    hold on
    if cfg==1
        legend(aa(fi), sprintf('F = %.2f MN', Fas(fi)*1e-6));
    end
end
aa1 = plot(UwxC{1}(end-2,:), (10.^UwxC{1}(end,:)).*sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UwxC{1}(1:end-3, :)).^2)), ...
    'k--', 'LineWidth', 1.5); hold on
aa2 = plot(Omegas(:, 1), sqrt(0.5)*abs(GM.NLTs.L*Phi{1}').*Qs', '-', 'Color', [1 1 1]*0.4,...
    'LineWidth', 1.5);
if cfg==1
    legend(aa1, 'EPMC')
    legend(aa2, 'RQNM')
end

set(gca, 'yscale', 'linear')
xlim([Wst Wen])
ylim([0 2.9])
xlabel('Frequency (rad/s)')
ylabel('RMS Amplitude (m)')

% inset
if cfg==1
    IS1 = [8200 8700 0 0.35];
elseif cfg==2
    IS1 = [8150 8400 0 0.75];
end

plot([IS1(1), IS1(1), IS1(2), IS1(2), IS1(1)], [IS1(3), IS1(4), IS1(4), IS1(3), IS1(3)], 'k-')
if cfg==1
    legend([aa(1:end); aa1; aa2], 'Location', 'northwest')
end

% is1 = axes('position', [0.18, 0.5, 0.25, 0.25]);
is1 = axes('position', [0.6, 0.425, 0.25, 0.25]);
for fi=1:length(Fas)
    plot(UCs{fi}(end,:), ...
        sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UCs{fi}(1:end-1, :)).^2)), ...
        '-', 'LineWidth', 2);
    hold on
end
plot(UwxC{1}(end-2,:), (10.^UwxC{1}(end,:)).*sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UwxC{1}(1:end-3, :)).^2)), ...
    'k--', 'LineWidth', 1.5); hold on
plot(Omegas(:, 1), sqrt(0.5)*abs(GM.NLTs.L*Phi{1}').*Qs', '-', 'Color', [1 1 1]*0.4,...
    'LineWidth', 1.5);

xlim([IS1(1), IS1(2)]);
ylim([IS1(3), IS1(4)]);
if cfg==1
    set(gca, 'XTick', [8250 8650])
    set(gca, 'YTick', [0.1 0.3])
elseif cfg==2
    set(gca, 'XTick', [8200 8350])
    set(gca, 'YTick', [0.1 0.6])
end

% print(sprintf('./FIGS/M1_FRESPWBB_cfg%d.eps', cfg), '-depsc')

%% Quantities for Synthesis
FEX = Fl(GM.Ndofs+(1:GM.Ndofs))+1j*Fl(2*GM.Ndofs+(1:GM.Ndofs));

Omsq_1 = @(famp) sqrt((Omegas(:, 1).^2-Cs(:,1).^2/2) + sqrt(Cs(:,1).^4/4-(Omegas(:,1).*Cs(:,1)).^2 + ((Phi{1}*FEX*famp)./Qs).^2));
Omsq_2 = @(famp) sqrt((Omegas(:, 1).^2-Cs(:,1).^2/2) - sqrt(Cs(:,1).^4/4-(Omegas(:,1).*Cs(:,1)).^2 + ((Phi{1}*FEX*famp)./Qs).^2));
AmpSc = sqrt(0.5)*abs(GM.NLTs.L*Phi{1}');
Phs = rad2deg(angle(Phi{1}*FEX));

%% Synthesis of FRF
figure(8)
clf()

aa = gobjects(length(Fas), 1);
for fi=1:length(Fas)
    aa1 = plot(UCs{fi}(end,:), ...
        sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UCs{fi}(1:end-1, :)).^2)), ...
        'b-', 'LineWidth', 2);
    hold on
%     if cfg==1
%         legend(aa(fi), sprintf('F = %.2f MN', Fas(fi)*1e-6));
%     end

    om1 = Omsq_1(Fas(fi));  
    om2 = Omsq_2(Fas(fi));
    
    ir1 = find(imag(om1)==0); 
    ir2 = find(imag(om2)==0);
    ir2 = ir2(end:-1:1);
    
    aa2 = plot([om1(ir1); om2(ir2)], [AmpSc(ir1).*Qs(ir1)' AmpSc(ir2).*Qs(ir2)'], 'r-', 'LineWidth', 1);
end
if cfg==1
    legend([aa1 aa2], 'HBM Reference', 'RQNM Synthesis', 'Location', 'northeast')
end

xlim([Wst Wen])

annotation('arrow', ([8400 8500]-Wst)/(Wen-Wst), [0.5 1]/3)
text(8400, 0.8, 'Increasing $||f_{ex}||$', 'fontsize', 16, 'interpreter', 'latex')

xlabel('Frequency (rad/s)')
ylabel('RMS Amplitude (m)')

print(sprintf('./FIGS/FRESPSYNTH_cfg%d.eps', cfg), '-depsc')

%% Phase
figure(9)
clf()

aa = gobjects(length(Fas), 1);
for fi=1:length(Fas)
    aa1 = plot(UCs{fi}(end,:), ...
        -rad2deg(angle(GM.NLTs.L*(UCs{fi}(GM.Ndofs+(1:GM.Ndofs),:)+1j*UCs{fi}(2*GM.Ndofs+(1:GM.Ndofs),:)))), ...
        'b-', 'LineWidth', 2);
    hold on
%     if cfg==1
%         legend(aa(fi), sprintf('F = %.2f MN', Fas(fi)*1e-6));
%     end

    om1 = Omsq_1(Fas(fi));  
    om2 = Omsq_2(Fas(fi));
    
    ir1 = find(imag(om1)==0); 
    ir2 = find(imag(om2)==0);
    ir2 = ir2(end:-1:1);
    
    aa2 = plot([om1(ir1); om2(ir2)], ...
        [Phs(ir1)-atan2d(Cs(ir1,1).*om1(ir1), (Omegas(ir1,1).^2-om1(ir1).^2)); Phs(ir2)-atan2d(Cs(ir2,1).*om2(ir2), (Omegas(ir2,1).^2-om2(ir2).^2))], ...
        'r-');
    hold on
%     aa2 = plot([om1(ir1); om2(ir2)], [AmpSc(ir1).*Qs(ir1)' AmpSc(ir2).*Qs(ir2)'], 'r-', 'LineWidth', 1);
end
% if cfg==1
%     legend([aa1 aa2], 'HBM Reference', 'RQNM Synthesis', 'Location', 'northeast')
% end
% 
xlim([Wst Wen])
ylim([-180 0])
set(gca, 'YTick', -180:45:0)
 
% annotation('arrow', ([8400 8500]-Wst)/(Wen-Wst), [0.5 1]/3)
% text(8400, 0.8, 'Increasing $||f_{ex}||$', 'fontsize', 16, 'interpreter', 'latex')
% 
xlabel('Frequency (rad/s)')
ylabel('Response Phase (degs)')

%% Plot 2-mode estimated parameters
load(sprintf('./DATA/Jenkins_RQNM_cfg%d.mat', cfg), 'As', 'Cs', 'Omegas', 'Phi', 'GM', 'Qs', 'UlC', 'Knl', 'Cnl', 'Q1s', 'Q2s');

figure(10)
clf()
set(gcf, 'outerposition', [100 20 1500 500])

for si=1:3
    subplot(1, 3, si)
%     set(gca, 'position', [0.13 0.11 0.2118 0.6364])

    imagesc(Qs, Qs, squeeze(Knl(si, :, :))')
    set(gca, 'YDir', 'normal')
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    xx = colorbar('northoutside');
    
%     set(gca, 'position', [0.13 0.01 0.2118 0.6364])
    
    switch si
        case 1
            xlabel(xx, '$k_{11}$', 'interpreter', 'latex', 'fontsize', 20)
            caxis([6.5 7.3]*1e7)
        case 2
            xlabel(xx, '$k_{12}$', 'interpreter', 'latex', 'fontsize', 20)
%             caxis([-5.7 1.5]*1e6)
            caxis([-5.7 5.7]*1e6)
        case 3
            xlabel(xx, '$k_{22}$', 'interpreter', 'latex', 'fontsize', 20)
            caxis([6 6.5]*1e8)
    end
    
    set(gca, 'XTick', [1e-2 1e0 1e2])
    set(gca, 'YTick', [1e-2 1e0 1e2])
    
    xlabel('$Q_1$')
    ylabel('$Q_2$')
    
%     set(gca, 'position', [0.13 0.2 0.2 0.5])
%     axis equal
end
print(sprintf('./FIGS/FitKvals_cfg%d.png', cfg), '-dpng')


figure(11)
clf()
set(gcf, 'outerposition', [100 600 1500 500])

for si=1:3
    subplot(1, 3, si)
    imagesc(Qs, Qs, squeeze(Cnl(si, :, :))')
    set(gca, 'YDir', 'normal')
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    xx = colorbar('northoutside');
    
    switch si
        case 1
            xlabel(xx, '$c_{11}$', 'interpreter', 'latex', 'fontsize', 20)
            caxis([65 425])
        case 2
            xlabel(xx, '$c_{12}$', 'interpreter', 'latex', 'fontsize', 20)
%             caxis([-1.5 212])
            caxis([-212 212])
        case 3
            xlabel(xx, '$c_{22}$', 'interpreter', 'latex', 'fontsize', 20)
            caxis([-600 600])
    end
    
    set(gca, 'XTick', [1e-2 1e0 1e2])
    set(gca, 'YTick', [1e-2 1e0 1e2])
    
	xlabel('$Q_1$')
    ylabel('$Q_2$')
%     axis equal
end
print(sprintf('./FIGS/FitCvals_cfg%d.png', cfg), '-dpng')