clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/MENGSHI_PFF/')
addpath('../../../RESEARCH/PFF_Code/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

analyze = false;
plotout = false;

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

%% Design Impulse
bwimp = 100;
fex = @(t) sin(bwimp*t).^2.*(t<=pi/bwimp)*0;
u0 = zeros(1,2)+V(:,1)'*1e1 + V(:,2)'*1e-1;  % 2.5e0, 1e-1
% u0 = rand(1,2)*1e0;
ud0 = zeros(1,2);
%% Conduct HHTA Simulation
fsamp = 32;
Ncyc = 400;
Tmax = Ncyc*2*pi/sqrt(Wr(2));

opt = struct('Display', 'waitbar');
famp = 5e6;
[T, u, ud, udd] = MDL.HHTAMARCH(0, Tmax, 1/fsamp, u0, ud0, @(t) fex(t)*famp*[1;0], opt);

%% Filter, then Hilbert-based postprocessing
Nfrt = 5e3;
x = [zeros(Nfrt,1); udd(1,:)'];
t = (0:length(x)-1)'/fsamp;
xfilt = [bandpass(x(end:-1:1), (Wr(1)+[-1 1]*0.5)/2/pi, fsamp) bandpass(x(end:-1:1), (Wr(2)+[-1 1]*0.4)/2/pi, fsamp)];
xfilt = xfilt(end:-1:1,:);

% %%
% [frqs, xf] = FFTFUN(t(1:end-Nfrt), x(Nfrt+1:end));
[frqs, xf] = FFTFUN(t, x);
[frqsf, xff] = FFTFUN(t(1:end-Nfrt-Nfrt/2+1), xfilt(Nfrt+Nfrt/2:end,:));

figure(100)
clf()
subplot(2,1,1)
plot(t, x, 'k-'); hold on
plot(t, xfilt(:,1), 'b-', t, xfilt(:,2), 'r-')
aa = gobjects(3,1);
aa(1) = plot(nan, nan, 'k-');
aa(2) = plot(nan, nan, 'b-', 'LineWidth', 2);
aa(3) = plot(nan, nan, 'r-', 'LineWidth', 2);
legend(aa, 'DOF 2', 'Mode 1 filt.', 'Mode 2 filt.')
xlabel('Time (s)')
ylabel('DOF 2 Acc. ($m/s^2$)')
subplot(2,1,2)
loglog(frqs*2*pi, abs(xf), 'k-'); hold on
loglog(frqsf*2*pi, abs(xff(:,1)), 'b-', 'LineWidth', 2)
loglog(frqsf*2*pi, abs(xff(:,2)), 'r-', 'LineWidth', 2)
xlim([1e-1 1e1])
xlabel('Frequency (rad/s)')
ylabel('DOF 2 Acc. ($m/s^2$)')

if plotout
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/RDOWN.png', '-dpng')
end

%%
tscale = 1000;
tm = t(1:end-Nfrt-Nfrt/2+1)/tscale;
xm1 = xfilt(Nfrt+Nfrt/2:end,1);
xm2 = xfilt(Nfrt+Nfrt/2:end,2);

save('./DATA/rdown.mat', 'tm', 'xm1', 'xm2', 'tscale');

%%
[Freq1,Damp1,Amp1]=FreqDampAmpExtraction(xm1,tm);
[Freq2,Damp2,Amp2]=FreqDampAmpExtraction(xm2,tm);

Freq1 = Freq1/tscale;
Freq2 = Freq2/tscale; 

% Amp1 = Amp1./((2*pi*Freq1).^2);
% Amp2 = Amp2./((2*pi*Freq2).^2);

%% 
figure(2)
clf()
subplot(2,1,1)
semilogx(Amp1, Freq1*2*pi)
subplot(2,1,2)
semilogx(Amp2, Freq2*2*pi)

save('./DATA/rdownproc.mat', 'tm', 'xm1', 'xm2', 't', 'x', 'xf', 'xfilt', 'xff', 'Amp1', 'Freq1', 'Damp1', 'Amp2', 'Freq2', 'Damp2')