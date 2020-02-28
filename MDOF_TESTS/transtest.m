clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/QUASISTATIC')
addpath('../ROUTINES/TRANSIENT')
addpath('../ROUTINES/SOLVERS')

model = 'BRB';

load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'K', 'M', 'R', 'L', 'Fv', 'TFM');

%% Load Mesh
Nds = dlmread(sprintf('../MODELS/%s/Nodes.dat', model));
Quad = dlmread(sprintf('../MODELS/%s/Elements.dat', model));

MESH = MESH2D(Nds, 3, [], Quad, 2);

%% Set Contact Function
MESH = MESH.SETCFUN(@(u, z, ud, P) ELDRYFRICT(u, z, ud, P, 0));

Pars = [1e12; 1e12; 1e12; 0.25];
pA = repmat(eye(4), MESH.Ne*MESH.Nq^2, 1);

U0 = zeros(size(K, 1), 1);

Z0 = zeros(2, MESH.Ne*MESH.Nq^2);
Prestress = 12000;

%% Fully Stuck Initial Guess 
Kstuck = zeros(size(L, 1));
Kstuck(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(1);
Kstuck(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(2);
Kstuck(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(3);

Kstuck = K + L'*Kstuck*L;

U0 = Kstuck\(Fv*Prestress);

Z0 = zeros(2, MESH.Ne*MESH.Nq^2);
%% Nonlinear Prestress Simulation
opts = struct('reletol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, 'etol', ...
              1e-6, 'ITMAX', 100, 'Display', true);
[Ustat, ~, eflag, ~, J0] = NSOLVE(@(U) QS_RESFUN([U; 0], Z0, Pars, L, pA, ...
						 MESH, M, K, Fv*Prestress, Fv*0), U0, opts);
[~, Z, ~, ~, ~, Txyn_p] = MESH.CONTACTEVAL(Ustat, Z0, Ustat*0, Pars, pA, L);  % Evaluate Slider
%% Linearized Modal Analysis
[V, D] = eigs(J0, M, 10, 'SM');
[D,si] = sort(diag(D));
Ws = sqrt(D)/2/pi;
V = V(:,si);
V = V./sqrt(diag(V'*M*V))';
% Rayleigh coefficients for damping
zts = [0.002; 0.003];
ab = [1./(2*2*pi*Ws([1 3])) 2*pi*Ws([1 3])/2]\zts;
C = ab(1)*M+ab(2)*J0;

% Excitation
freq = 500;
fex = @(t) R(3, :)'*(400*sin(2*pi*freq*t)^2*(t<0.5/freq))+Fv*Prestress;

% Initial Conditions
U0 = Ustat;
Ud0 = Ustat*0;
Z0 = Z;
disp('INITIAL CONDITIONS SET');

%% HHTA Hysteretic
ABG = [0, 1/4, 1/2];  % Unconditionally Stable Newmark-Alpha
%% ABG = [0, 1/6, 1/2];  % Implicit linear acceleration
%% ABG = [-0.1, 1/6, 1/2];  % HHT-Alpha

T0 = 0;
T1 = 2.5;
dT = 1e-4;  % 2000 Hz Nyquist

opts = struct('reletol', 1e-12, 'etol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, 'Display', true, 'ITMAX', 100);
[Th, Xh, zh, Xdh, Xddh] = HHTA_NONLIN_HYST(M, C, K, fex,
					   @(t, x, z, xd) MESH.CONTACTEVAL(x, z, xd, Pars, pA, L),
					   U0, Z0, Ud0, T0, T1, dT, ABG(1), ABG(2), ABG(3), opts);

plot(Th, R(3, :)*Xh, 'o-', 'LineWidth', 2)

save('./DATS/RUN1.mat', 'Th', 'Xh', 'zh', 'Xdh', 'Xddh', '-v7')

% [Thh, Xhh, zhh, Xdhh, Xddhh] = HHTA_NONLIN_HYST(M, C, K, fex,
% 						@(t, x, z, xd) MESH.CONTACTEVAL(x, z, xd, Pars, pA, L),
% 						Xh(:,end), zh(:,:,end), Xdh(:,end), T1, 2*T1, dT, ABG(1), ABG(2), ABG(3), opts);
% plot(Thh, R(3, :)*Xhh, 'o-')
