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
C = 0.01*K;  % 0.01*K

MDL = MDOFGEN(M, K, C, eye(2));

kt = 4;  % 4
muN = 2;  % 2
MDL = MDL.SETNLFUN(2+3, [1 0], @(t, u, varargin) JENKNL(t, u, kt, muN, varargin{:}), [], 4);
K0 = K+kt*[1 0;0 0]; % Linearized stiffness
[V, Wr] = eig(K0, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

disp(diag(V'*C*V)./(2*Wr))

%% Design Impulse
bwimp = 100;
fex = @(t) sin(bwimp*t).^2.*(t<=pi/bwimp)*1e-1*0;
u0 = zeros(1,2)+V(:,1)'*2.5 + V(:,2)'*1;  % 2.5e0, 1e-1
% u0 = [1;0]*2;
% u0 = rand(1,2)*0e0;
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
tscale = 1000;
t = t/tscale;

save('./DATA/rdowndat_eldrfr.mat', 't', 'x', 'tscale')

%% 
[frqs, xf] = FFTFUN(t, x);

figure(1)
clf()
subplot(2,1,1)
% plot(t*tscale, x)
plot(T, u(1,:))
xlabel('Time (s)')
zlabel('DOF 1 Acceleration')
subplot(2,1,2)
semilogy(2*pi*frqs/tscale, abs(xf))
xlabel('Frequency (Hz)')
xlim([0 5])
zlabel('DOF 2 Acceleration')