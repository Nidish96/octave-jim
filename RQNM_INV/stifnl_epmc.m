clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

%% Setup System
M = [1 0;0 1];
K = [2 -1;-1 2];
L = [1 0];
% Lf = [0.5; 0];
% fnl = @(t, u, ud) deal(u.^3, 3*u.^2, zeros(size(u)));
fnl = @(t, u, ud) deal(0.5*u.^3, 1.5*u.^2, zeros(size(u)));
%% Setup Model
GM = MDOFGEN(M, K, M*0, eye(2));
% GM = GM.SETNLFUN(1+5, L, fnl, Lf);
GM = GM.SETNLFUN(1+3, L, fnl);

%% Linear Modes
[V, Wsp] = eig(K, M);
[Wsp,si] = sort(sqrt(diag(Wsp)));
V = V(:,si);
V = V./sqrt(diag(V'*M*V))';

%% EPMC
h = 0:7;
Nhc = sum((h==0)+2*(h~=0));
Nt = 2^10;

As = -2;
Ae = 3;
da = 0.01;

Copt = struct('Nmax', 5000, ...
    'dsmax', 0.1, 'dsmin', da/1000, ...
    'angopt', 1e-2);

Uwx0 = [kron([0; 1; 1; zeros(Nhc-3,1)], V(:,1)); Wsp(1); 0];
Fl = kron([0; 1; 0; zeros(Nhc-3,1)], L');

[UwxC, dUwxC] = CONTINUE(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt), Uwx0, As, Ae, da, Copt);

%% Save
save('./DATA/Stifnl_EPMC.mat', 'UwxC', 'dUwxC', 'h', 'Nhc');

%% Load
load('./DATA/Stifnl_EPMC.mat', 'UwxC');

%% Plot
figure(1);
clf();
plot(10.^UwxC(end,:), UwxC(end-2,:), '.-'); hold on

set(gca, 'XScale', 'log')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (rad/s)')

figure(2);
clf()
plot(10.^UwxC(end,:), UwxC(end-1,:)./(2*UwxC(end-2,:)), '.-'); hold on

set(gca, 'XScale', 'log')
xlabel('Modal Amplitude')
ylabel('Damping Factor')