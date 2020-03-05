clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/QUASISTATIC')
addpath('../ROUTINES/TRANSIENT')
addpath('../ROUTINES/HARMONIC')
addpath('../ROUTINES/SOLVERS')

model = 'BRB';

load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'K', 'M', 'R', ...
     'L', 'Fv', 'TFM');
prec = 'single';
K = single(K);
M = single(M);
R = single(R);
Fv = single(Fv);
clear TFM

%% Load Mesh
Nds = single(dlmread(sprintf('../MODELS/%s/Nodes.dat', model)));
Quad = single(dlmread(sprintf('../MODELS/%s/Elements.dat', model)));
Ne = size(Quad,1);

Nq = 1;
MESH = MESH2D(Nds, 3, [], Quad, Nq);
MESH = MESH.SETCFUN(@(u, z, ud, P) ELDRYFRICT(u, z, ud, P, 0), zeros(2, Ne*Nq^2));  % Contact Function for static analysis
MESH = MESH.SINGLEPREC();

% Parameterization
Pars = [1e12; 1e12; 1e12; 0.25];
pA = sparse(repmat(eye(4), MESH.Ne*MESH.Nq^2, 1));

%% PRESTRESS ANALYSIS
Prestress = single(12000);

% Fully Stuck Initial Guess 
Kstuck = sparse(size(L, 1), size(L, 1));
Kstuck(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(1);
Kstuck(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(2);
Kstuck(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(3);

Kstuck = K + L'*Kstuck*L;

U0 = Kstuck\(Fv*Prestress);
Z0 = zeros(2, MESH.Ne*MESH.Nq^2, 'single');
%% Nonlinear Prestress Simulation
opts = struct('reletol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, 'etol', ...
              1e-6, 'ITMAX', 100, 'Display', true, 'Dscale', ones(size(U0), 'single'));
% opts.Dscale = ones(size(U0))*max(abs(U0));

[Ustat, ~, eflag, ~, J0] = NSOLVE(@(U) QS_RESFUN([U; 0], MESH.z, Pars, L, pA, ...
						 MESH, M, K, Fv*Prestress, Fv*0), U0, opts);
[~, Z, ~, ~, ~, Txyn_p] = MESH.CONTACTEVAL(Ustat, Z0, Ustat*0, Pars, pA, L);  % Evaluate Slider
disp(sum(MESH.Tm*Txyn_p(3,:)')/(Prestress*3))

%% Linearized Modal Analysis
[V, D] = eigs(double(J0), double(M), 10, 'SM');
[D,si] = sort(diag(D));
Ws = sqrt(D)/2/pi;
V = V(:,si);
V = V./sqrt(diag(V'*M*V))';

% SET LINEAR DAMPING
zts = [0.002; 0.003];
ab = [1./(2*2*pi*Ws([1 3])) 2*pi*Ws([1 3])/2]\zts;
C = ab(1)*M+ab(2)*J0;

%% HARMONIC BALANCE
Nt = uint32(128);
h = uint32([0 1]);  Nhc = uint32(sum(h==0)+2*sum(h~=0));
Nd = uint32(size(K, 1));

% Linear Forcing
fa = single(5);
Fl = single(kron([0 fa 0 zeros(1,Nhc-3,'single')], R(3,:))');
if h(1)~=0
  Fl(1:Nd) = [];
else
  Fl(1:Nd) = Fv*Prestress;
end

wfrc = single(2*pi*150);

Elin = HARMONICSTIFFNESS(M, C, J0, wfrc, h(h~=0));
if h(1)~=0
  U0 = single(Elin\double(Fl));
else
  U0 = single([Ustat; Elin\double(Fl(Nd+1:end))]);
end

opts = struct('reletol', 1e-10, 'rtol', 1e-6, 'utol', 1e-6, 'etol', ...
              1e-6, 'ITMAX', 100, 'Display', true, 'Dscale', ones(size(U0), 'single'));

MESH = MESH.SETCFUN(@(U, h, Nt, P) ELDRYFRICT_HB(U, h, Nt, P, ...
                                                 single(0)), []);  % Contact Function

tic
[R0, dR0] = MDOF3D_NLHYST_HBRESFUN([U0; wfrc], Pars, L, pA, MESH, ...
                                   M, C, K, Fl, h, Nt, 1:MESH.Nn*MESH.dpn);
toc

if h(1)==0
  opts.Dscale(1:Nd) = ones(Nd,1)*max(abs(Ustat));
  opts.Dscale(Nd+1:end) = ones(Nd*(Nhc-1),1)*max(abs(U0(Nd+1:end)));

  U0 = U0./opts.Dscale;
end

[Uws, ~, eflag] = NSOLVE(@(U) MDOF3D_NLHYST_HBRESFUN([U; wfrc], Pars, L, pA, MESH, M, C, K, Fl, h, Nt, 1:MESH.Nn*MESH.dpn), U0, opts);

				% Sweep
Wsweep = 2*pi*[140:2:155 155:0.5:162 163:2:180];
Nsweep = length(Wsweep);

Xsweep = zeros(Nhc, Nsweep);
Xsweep(:, 1) = (R(3,:)*reshape(Uws, Nd, Nhc))';

opts.ITMAX = 20;
for iw=iw:Nsweep
  [Uws, ~, eflag] = NSOLVE(@(U) MDOF3D_NLHYST_HBRESFUN([U; Wsweep(iw)], Pars, L, pA, MESH, ...
						     M, C, K, Fl, h, Nt, 1:MESH.Nn*MESH.dpn), ...
			   Uws, opts);
  Xsweep(:, iw) = (R(3,:)*reshape(Uws, Nd, Nhc))';
  fprintf('Done %d/%d W=%f\n', iw, Nsweeps, Wsweep(iw)/2/pi);
end

save('./DATS/RUN1_HB.mat', 'Xsweep', 'Wsweep', 'fa');

% 				% CONTINUATION
% Copt = struct('Nmax', 50, 'Display', 1, 'angopt', 1e-2);
% Ubws = CONTINUE(@(Uw) MDOF3D_NLHYST_HBRESFUN(Uw, Pars, L, pA, MESH, M, C, K, Fl, h, Nt, 1:MESH.Nn*MESH.dpn), Uws, 2*pi*150, 2*pi*170, 2*pi*2, Copt);

