clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/QUASISTATIC')
addpath('../ROUTINES/TRANSIENT')
addpath('../ROUTINES/SOLVERS')
addpath('../ROUTINES/HARMONIC')
addpath('../ROUTINES/MENGSHI_PFF')

model = 'BRB';
iR = 3;

runs = {'1_Low', '2_Med', '3_High'};

load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'R');
%% Time Domain Plots

aa = gobjects(3,1);
bb = gobjects(3,1);
for iR=1:3
    load(sprintf('./DATS/RUN%s.mat',runs{iR}), 'Th', 'Xddh', 'Finput', ...
         'famp')
    figure(1)
    if iR==1
        clf()
    end
    
    aa(iR) = plot(Th, R(3,:)*Xddh/famp); hold on
    legend(aa(iR), sprintf('F = %d N', famp));
    
    figure(2)
    if iR==1
        clf()
    end
    [freqs, Xf] = FFTFUN(Th', (R(3,:)*Xddh)');
    bb(iR) = plot(freqs, abs(Xf)/famp, 'LineWidth', 2); hold on
    legend(aa(iR), sprintf('F = %d N', famp));
end
figure(1)
legend(aa(1:3), 'Location', 'northeast')

xlim([0 1])
xlabel('Time (t)')
ylabel('Acceleration/Force (m s^{-2}/N)')

print('./FIGS/RINGDOWN.eps', '-depsc')

figure(2)
legend(bb(1:3), 'Location', 'northeast')
set(gca, 'Yscale', 'log')
xlim([100 200])
% ylim([1e-4 1e-1])
xlabel('Frequency (Hz)')
ylabel('FRF Amplitude (m s^{-2}/N)')

print('./FIGS/TFRFM1.eps', '-depsc')

%% Backbone plot
load('./DATS/TRANSPROC.mat', 'Amp', 'Freq', 'Damp')

Freqs = smoothdata(Freq);
Damps = smoothdata(Damp);

ic = find((Amp(1:end-1)-2.5e-1).*(Amp(2:end)-2.5e-1)<0);
Amps = Amp(1:ic);
Freqs = Freqs(1:ic);
Damps = Damps(1:ic);

save('./DATS/TRANSPROC.mat', 'Amp', 'Freq', 'Damp', 'Amps', 'Freqs', 'Damps')

figure(3)
clf()
subplot(2, 1,1)
semilogx(Amp, Freq, 'k.'); hold on
semilogx(Amps, Freqs, 'b-', 'LineWidth', 2); hold on
ylim([154 161])

ylabel('Natural Frequency (Hz)')

subplot(2, 1,2)
semilogx(Amp, Damp, 'k.'); hold on
semilogx(Amps, Damps, 'b-', 'LineWidth', 2); hold on

ylim([0 0.025])
ylabel('Damping Factor')
xlabel('Acceleration Amplitude (m s^{-2})')

print('./FIGS/BACKBONE.eps', '-depsc')