clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/QUASISTATIC')
addpath('../ROUTINES/TRANSIENT')
addpath('../ROUTINES/SOLVERS')
addpath('../ROUTINES/MENGSHI_PFF')

model = 'BRB';
iR = 3;

runs = {'1_Low', '2_Med', '3_High'};

load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'R');
load(sprintf('./DATS/RUN%s.mat',runs{iR}), 'Th', 'Xddh', 'Finput', ...
     'famp')

md = 1;  % Mode to analyze

switch md
  case 1
    fn = 159;  % central frequency in filtering (Hz)
    bw = 90; % bandwidth of frequency in the bandpass filter
    n_bp = 4; % order of the bandpass filter
    sca = 1;
end

%% Input (need to be changed properly)
TimeDuration    = Th(end); % the time length to be analyzed 
fs 	        = 1.0/(Th(2)-Th(1));
% NOTE: TimeDuration depends on noise (not too long, can be changed
% according to the calculated results below)

%% Raw Data Filtering 

data_raw = R(3, :)*Xddh;

[xR, ~, t] = RevForFilt(data_raw, Th, fn, n_bp, bw);
[ X, X_abs, Fz1] = FFT_simp(data_raw, fs);
[ Y, Y_abs, Fz2] = FFT_simp(xR, fs);

figure(1) % FFT of raw data and filtered data
clf()
semilogy(Fz1,X_abs)
hold on;
plot(Fz2,Y_abs,'r')
xlim([0,3000])

% xlim([0,250])
% ylim([1e-1 1e1])
legend('Raw Signal','Filtered Signal', 'Location', 'northwest')
xlabel('Frequency (Hz)', 'interpreter', 'latex')
ylabel('Acceleration Amplitude (m s$^{-2}$)', 'interpreter', ...
       'latex')
%% Cut Data Length
temp1=find(Th<=TimeDuration);temp1=temp1(end);
xR=xR(1:temp1);t=Th(1:temp1);

%% Calculation
Ln = length(t);
LnA = fix(Ln/sca);
[Freq,Damp,Amp,Time] = FreqDampAmpExtraction(xR(1:LnA),t(1:LnA));
%% Figure for Frequency (Backbone Curve)
figure(3)
if iR==1
    clf()
end
subplot(2, 1, 1)
semilogx(Amp./(2*pi*Freq).^2/9.81,Freq,'.'); hold on
% xlabel('Amplitude (m)');
ylabel('Frequency (Hz)', 'interpreter', 'latex');

xlim([3e-8 1e-4])

%% Figure for Damping
subplot(2, 1, 2)
semilogx(Amp./(2*pi*Freq).^2/9.81,Damp*100,'.'); hold on
set(gca, 'YScale', 'log')
xlabel('Amplitude (m)', 'interpreter', 'latex');
ylabel('Damping Ratio (%)', 'interpreter', 'latex');
ylim([0 max(Damp*120)])

xlim([3e-8 1e-4])

%% Save data
if iR==3
    save('./DATS/TRANSPROC.mat', 'Amp', 'Freq', 'Damp');
end
% figure(1)
% if iR==1
%     clf()
% end
% plot(Th, R(3,:)*Xddh/famp); hold on
