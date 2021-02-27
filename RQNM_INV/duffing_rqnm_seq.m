clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/FEM/')

%% Setup system (1D Bar)
ell = 1.0;
A = pi*1e-4;
E = 210e9;
rho = 7800;
Ne = 8;

Le = ell/Ne;
Me = 1/3*[2 1;1 2]*Le/2*rho*A;
Ke = 1/2*[1 -1;-1 1]*2/Le*E*A;

M = zeros(Ne+1);
K = zeros(Ne+1);
for e=1:Ne
    M(e:e+1,e:e+1) = M(e:e+1,e:e+1) + Me;
    K(e:e+1,e:e+1) = K(e:e+1,e:e+1) + Ke;
end
Lb = eye(Ne+1);
Lb(:, 1) = [];  % Boundary Condition (DOF 1 fixed)

Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;

%% Damping Factors
Zetas = [0.4; 0.4]*1e-2;

[V, Wsp] = eig(Kb, Mb);
[Wsp, si] = sort(sqrt(diag(Wsp)));
V = V(:, si);
V = V./diag(sqrt(V'*Mb*V))';

if Ne>2
    ab = [1./(2*Wsp(1:length(Zetas))) Wsp(1:length(Zetas))/2]\Zetas;
    Cb = ab(1)*Mb + ab(2)*Kb;
else
    Cb = 2*Zetas(1)/Wsp(1)*Kb;
end

%% Setup model
GM = MDOFGEN(Mb, Kb, Cb, Lb);

kc = 1e6;
fnl = @(t,u,ud) deal(kc*u.^3, 3*kc*u.^2,zeros(size(u)));
GM = GM.SETNLFUN(1+3, Lb(end,:), fnl);

%%
Q = 1e1;

opt = optimset('Jacobian', 'on', 'Display', 'iter');
Uls = fsolve(@(ul) GM.RQMRESFUN([ul; Q],0), [V(:, 1)*Q; Wsp(1)^2], opt)

%% Sequential Implementation
Q = 1e5;

u0 = V(:,1)*Q;
l0 = Wsp(1)^2;

ITMAX = 100;
CSs = zeros(ITMAX, 1);

opt.Display = 'off';
for it=1:ITMAX
    % 1. Update u
    % Linear Newton Update
    [Rs, Js] = GM.STATRESFUN(u0, l0*GM.M*u0);

    u1 = u0 - Js\Rs;
    du1dlam = Js\(GM.M*u0);

    % Nonlinear Solve Update
%     u1 = fsolve(@(u) GM.STATRESFUN(u, l0*GM.M*u0), u0, opt);
%     [Rs, Js] = GM.STATRESFUN(u1, l0*GM.M*u0);
%     du1dlam = Js\(GM.M*u0);

    % 2. Update lam
    c = u1'*GM.M*u1/Q^2-1;
    dcdlam = 2*u1'*GM.M*du1dlam/Q^2;

    l1 = l0 - c/dcdlam;
    
    % Update variables
    l0 = l1;
    u0 = u1*Q/sqrt(u1'*GM.M*u1);
    
    % Accumulator & print
    CSs(it) = c;
    fprintf('%d %e\n', it, c)
end

figure(1);
% clf()
hold on
semilogy(1:ITMAX, abs(CSs), 'x-')