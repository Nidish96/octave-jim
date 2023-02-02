clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

analyze = false;
plotout = false;

%% Parameters
M = eye(2);
K = [1 -1;-1 2];
W = sqrt(eig(K,M));

kt = 1;  % 4
muN = 1e-1;  % 2
K0 = K+kt*[1 0;0 0]; % Linearized stiffness
% C = 5e-3*K0;  % 0.01*K
C = 1e-2*M;

MDL = MDOFGEN(M, K, C, eye(2));
MDL = MDL.SETNLFUN(2+3, [1 0], @(t, u, varargin) JENKNL(t, u, kt, muN, varargin{:}), [], 4);
[V, Wr] = eig(K0, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

zts = diag(V'*C*V)./(2*Wr);
disp(zts)
disp([W Wr])

%% Design Impulse
bwimp = 5;
fex = @(t) sin(bwimp*t).^2.*(t<=pi/bwimp);
% u0 = zeros(1,2)+V(:,1)'*2.5 + V(:,2)'*1;  % 2.5e0, 1e-1
% u0 = [0;0]*2;
% u0 = rand(1,2)*0e0;
u0 = zeros(1,2);
ud0 = zeros(1,2);
%% Conduct HHTA Simulation
fsamp = 16;
Ncyc = 400;
Tmax = Ncyc*2*pi/Wr(2);

opt = struct('Display', 'waitbar');
famp = 1e1;
[T, u, ud, udd] = MDL.HHTAMARCH(0, Tmax, 1/fsamp, u0, ud0, @(t) fex(t)*famp*[1;0], opt);

%% Filter, then Hilbert-based postprocessing
Nfrt = 5e3;
x = [zeros(Nfrt,1); u(1,:)'];
xd = [zeros(Nfrt,1); ud(1,:)'];
xdd = [zeros(Nfrt,1); udd(1,:)'];
t = (0:length(x)-1)'/fsamp;
tscale = 1e3;
t = t/tscale;

save('./DATA/rdowndat_eldrfr.mat', 't', 'x', 'xd', 'xdd', 'tscale')

%% 
[frqs, xf] = FFTFUN(t, [x xd xdd]);
[~, ff] = FFTFUN(t, fex(t*tscale));

figure(1)
clf()
stackedplot(t*tscale, [x xd xdd [zeros(Nfrt,1); fex(T(:))]], '.-')

figure(2)
clf()
h = stackedplot(2*pi*frqs/tscale, abs([xf ff]));
% for i=1:4; h.AxesProperties(i).YScale = 'log'; end
% xlim([0 bwimp])

figure(3)
clf()
semilogy(2*pi*frqs/tscale, abs(xf(:,3)./ff(:,1)))
yyaxis right
semilogy(2*pi*frqs/tscale, abs(ff))
xlim([0 4*bwimp])