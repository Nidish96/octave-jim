clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/QUASISTATIC')
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

[F, Z, ~, ~, ~, Txyn_p] = MESH.CONTACTEVAL(L*Ustat, Z0, Ustat*0, Pars, pA);  % Evaluate Nodal Forces

figure(1)
clf()
MESH.SHOWFIELD2D(Txyn_p(3,:)')

%% Linearized Modal Analysis
[V, D] = eigs(J0, M, 10, 'SM');
[D,si] = sort(diag(D));
Ws = sqrt(D)/2/pi;
V = V(:,si);
V = V./sqrt(diag(V'*M*V))';

%% RQNM
opts = struct('reletol', 1e-10, 'rtol', 1e-6, 'utol', 1e-6, 'etol', ...
              1e-6, 'ITMAX', 100, 'Display', true);

q = 1e-3;
U0 = [Ustat+q*V(:,1);D(1)];

% U0 = Urq;
[Urq, ~, eflag] = NSOLVE(@(Ul) RQNMA_RESFUN([Ul; q], Z0, Pars, L, pA, MESH, M, K, Fv*Prestress, ...
                                            Ustat), U0, opts);
disp([Ws(1) sqrt(Urq(end))/2/pi])
