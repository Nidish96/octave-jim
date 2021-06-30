clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

% TAmps = {'Low', 'Med', 'High'};
TAmps = {'Low', 'High'};
Modes = [1 2 3];

mi = 2;
%%
cols = DISTINGUISHABLE_COLORS(length(TAmps));
aa = gobjects(size(TAmps));

figure((mi-1)*2+1)
clf()
figure((mi-1)*2+2)
clf()
for ai=1:length(TAmps)
    load(sprintf('MATFILES/Mode%d_%s.mat', Modes(mi), TAmps{ai}), 'AMP_avg', 'DAM_avg', 'FRE_avg', ...
        'Amps', 'Freqs', 'Damps');
    % Gaussian Process Regression
%     if mi==3
%         mamps = cell2mat(Amps')./(2*pi*cell2mat(Freqs')).^2; 
%         mfreqs = cell2mat(Freqs');
%         mdamps = cell2mat(Damps');
%         AMP_avg = logspace(log10(min(mamps)), log10(max(mamps)), 1000)';
%         mdlW = fitrgp(log10(mamps), mfreqs);
%         mdlZ = fitrgp(log10(mamps), log10(abs(mdamps)));
%         FRE_avg = predict(mdlW, log10(AMP_avg));
%         DAM_avg = 10.^predict(mdlZ, log10(AMP_avg));
%         
%         AMP_avg = AMP_avg.*(2*pi*FRE_avg).^2;
%     end
    
    figure((mi-1)*2+1)
    for i=1:length(Amps)
        semilogx(Amps{i}./(2*pi*Freqs{i}).^2, Freqs{i}, '-', 'Color', 0.6*[1 1 1]); hold on
    end
    aa(ai) = semilogx(AMP_avg./(2*pi*FRE_avg).^2, FRE_avg, '-', 'LineWidth', 2, 'Color', cols(ai,:));
%     aa(ai) = semilogx(mx, mW, '-', 'LineWidth', 2, 'Color', cols(ai,:));
    legend(aa(ai), sprintf('Test Amplitude: %s', TAmps{ai}));
    
    figure((mi-1)*2+2)
    for i=1:length(Amps)
        semilogx(Amps{i}./(2*pi*Freqs{i}).^2, Damps{i}*100, '-', 'Color', 0.6*[1 1 1]); hold on
    end
    semilogx(AMP_avg./(2*pi*FRE_avg).^2, DAM_avg*100, '-', 'LineWidth', 2, 'Color', cols(ai,:));
%     semilogx(mx, mZ*100, '-', 'LineWidth', 2, 'Color', cols(ai,:));
end
figure((mi-1)*2+1)
set(gcf, 'Color', 'white')
ll=legend(aa, 'Location', 'northeast');
% legend(aa, 'Location', 'southwest')
if mi~=3
    set(ll, 'Visible', 'off')
end
xlabel('Response Amplitude (m)')
ylabel('Natural Frequency (Hz)')
xlim([1e-9 2e-5])
export_fig(sprintf('./FIGS/Mode%d_ExpbbW.eps', Modes(mi)), '-depsc')

figure((mi-1)*2+2)
set(gcf, 'Color', 'white')
xlabel('Response Amplitude (m)')
ylabel('Damping Factor (\%)')
xlim([1e-9 2e-5])
export_fig(sprintf('./FIGS/Mode%d_ExpbbZ.eps', Modes(mi)), '-depsc')