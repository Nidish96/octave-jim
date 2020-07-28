clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

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

kc = 1e7;
fnl = @(t,u,ud) deal(kc*u.^2, 2*kc*u, zeros(size(u)));
GM = GM.SETNLFUN(1+3, Lb(end,:), fnl);

%% Check
h = 0:5;
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^9;

Fl = kron([0; 1; 0; zeros(Nhc-3,1)], Lb(end,:)');

As = -2;
Ae = 1.5;
da = 0.01;

Uwx0 = [kron([0; 0; 1; zeros(Nhc-3,1)], V(:,1)); Wsp(1); 2*Zetas(1)*Wsp(1)];
Dscale = [1e-2*ones(Nhc*GM.Ndofs,1); Wsp(1); 2*Zetas(1)*Wsp(1); 1.0];

Copt = struct('Nmax', 500, 'Dscale', Dscale, 'dsmax', 0.05);
[UwxC, dUwxC] = CONTINUE(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt), Uwx0, As, Ae, da, Copt);

%% Save
save('./DATA/Quadrat_EPMC.mat', 'UwxC', 'dUwxC', 'h', 'Nhc');

%% Load
load('./DATA/Quadrat_EPMC.mat', 'UwxC');

%% Plot
figure(1);
clf();
plot(10.^UwxC(end,:), UwxC(end-2,:)); hold on

set(gca, 'XScale', 'log')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (rad/s)')

figure(2);
clf()
plot(10.^UwxC(end,:), UwxC(end-1,:)./(2*UwxC(end-2,:))); hold on

set(gca, 'XScale', 'log')
xlabel('Modal Amplitude')
ylabel('Damping Factor')