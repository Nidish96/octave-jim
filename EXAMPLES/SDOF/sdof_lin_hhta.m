clc
clear all
addpath('../../ROUTINES/')

%%
par.m = 1.0;
par.c = 0.2;
par.k = (2*pi)^2;

x0 = [0;0];

MDL = MDOFGEN(par.m, par.k, par.c, 1.0);

%%
fsamp = 32;
blocksize = 2048;

finp = 1.0;

opt = struct('Display', 'waitbar', 'alpha', 0, 'beta', 1/6, 'gamma', 1/2);
[T, u, ud, udd] = MDL.HHTAMARCH(0, (blocksize-1)/fsamp, 1/fsamp, x0(1), x0(2), @(t) cos(2*pi*finp*t), opt);

%% Lsim
A = [0 1.0;
    -par.k/par.m -par.c/par.m];
B = [0; 1/par.m];
C = [1 0];
D = 0;
sys = ss(A,B,C,D);

Y = lsim(sys, cos(2*pi*finp*T), T, x0);

%%
figure(1)
clf()
plot(T, u); hold on
plot(T, Y, '.')