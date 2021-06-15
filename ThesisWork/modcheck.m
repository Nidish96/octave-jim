clc
clear all

%% Params
M = eye(3);
K = (2*pi*diag([170 180 2])).^2;

cn = 9e9;
ct = 1.3e6;
mu = 0.2840;
lam = 9.4037e5;
gap = 7.0260e-4;

Fv = [1e-5; 1e-5; 1e-2];
Prestress = 12e3;

%% Contact Model
fnl = @(t, u, varargin) EXPROUGHFRICT(t, u, ct, cn, lam, mu, gap, varargin{:});

%% Create Object
GM = MDOFGEN(M, K, zeros(size(M)), eye(3));

% Kt = [1e12; 1e12; 0];
% kn = 1e12;
% mu = 0.25;
% gap = 0;
% 
% fnl = @(t, u, varargin) ELDRYFRICT(t, u, Kt, kn, mu, gap, varargin{:});
GM = GM.SETNLFUN(2+5, ...
    eye(3), ...
    fnl, ...
    eye(3));

%% Linearized Approximation
K0 = diag([ct/cn*Prestress*Fv(3);ct/cn*Prestress*Fv(3);lam*Prestress*Fv(3)]);

%% Prestress
opts = struct('reletol', 1e-6, 'Display', true, 'lsrch', 0, 'ITMAX', 100);

U0 = zeros(size(K,1),1);
U0 = (K+K0)\(Fv*Prestress);
% U0(3:3:MESH.Nn) = U0(3:3:MESH.Nn)+gapo;
% U0 = L(3:3:MESH.Nn*3, :)'*(MESH.Qm\(gap+log(Prestress./Aels./cn)./lam));

% U0(3) = gap;
U0(3) = (log(Prestress*Fv(3))-log(cn))/lam+gap;
[Ustat, ~, ~, ~, Jstat] = NSOLVE(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, opts);  %% FIX GRADIENT !!!
% U0 = Ustat;

%%
% U0(3) = gap;
% U0 = Ustat;
opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display', 'iter');
Ustat = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, opts);
%%
opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter', 'CheckGradients', true);
Ustat = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, opts);

%%
[R, J] = GM.STATRESFUN(Ustat, Fv*Prestress);

%% Linearization
[V, D] = eig(J, M);
[Ws, si] = sort(sqrt(diag(D)));
V = V(:, si);

%% RQNM
qamp = 1e12;

Ul0 = [Ustat+V(:,1)*qamp; Ws(1)^2];
opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter', 'CheckGradients', false);
Us = fsolve(@(Ul) GM.RQMRESFUN([Ul; qamp], 0, Fv*Prestress, Ustat), Ul0, opts)

% UlC(1:end-1, Na+ia) = fsolve(@(ul) GM.RQMRESFUN([ul; As(Na+ia)], 0, Fv*Prestress, Ustat), ul0, fopts);

%%
GM.RQMRESFUN([Us; qamp], 0, Fv*Prestress, Ustat)