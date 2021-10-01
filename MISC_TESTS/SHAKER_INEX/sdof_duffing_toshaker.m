clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/HARMONIC/')

%% SDOF Duffing Oscillator Model
m = 3.6;
k = (2*pi*150)^2*m;
c = 2*0.002*sqrt(k/m)*m;
b = 1e17;

MDL = MDOFGEN(m, k, c, 1.0);
fnl = @(t,u,ud) deal(b*u.^3, 3*b*u.^2, zeros(size(u)));
MDL = MDL.SETNLFUN(1+3, 1, fnl);

%% Shaker Mode (State-Space)
%%%%% mechanical parameters
M_T = 0.0243; % Table Masse [kg]
M_C = 0.0190; % Coil Masse [kg]

K_C = 8.4222*10^7; % spring constant Coil-Table [N/m]
C_C = 57.1692; % damping constant Coil-Table [Ns/m]

K_T = 20707; % spring constant Table-Ground [N/m]
C_T = 28.3258; % damping constant Table-Ground [Ns/m]

%%%%% electrical
L = 140*10^-6; % inductivity [H]
R = 3.00; % Resistance [Ohm]
Tau = 15.4791; % shaker ratio of thrust to coil current [N/A]

%%%%% State space model
E_shaker = diag([L, 1, M_C, 1, M_T]);
A_shaker = [- R 0 -Tau 0 0; 0 0 1 0 0; ...
    Tau -K_C -C_C K_C C_C; 0 0 0 0 1; ...
    0 K_C C_C -(K_C+K_T) -(C_C+C_T)];
B_shaker = [1 0 0 0 0; 0 0 0 0 1]';  % First column is voltage input, second for input noise
C_shaker = [0 0 0 1 0; 0 0 0 0 1];
D_shaker = [0 0 0 0];

% stinger parameters
E_stinger = 210000 ;%in N/mm^2
A_stinger = pi*2^2; % in mm
l_stinger = 0.0200; %in m
k_stinger = (E_stinger*A_stinger)/l_stinger;

% 2nd Order Shaker Model
M_SH = diag([0 M_C M_T]);
C_SH = [L Tau 0; 0 C_C -C_C; 0 -C_C C_C+C_T];
K_SH = [R 0 0;-Tau K_C -K_C;0 -K_C K_C+K_T];

Nshaker = [0; 0; 1];
Stinger_K = diag([0 0 k_stinger]);

% [MDL, Fshape] = MDL.ATTACHSHAKER(E_shaker, A_shaker, B_shaker, Shaker_K, Shaker_K*0, Nshaker);
MDL = MDL.ATTACHSHAKER2(M_SH, C_SH, K_SH, Stinger_K, Stinger_K*0, Nshaker);
Fshape = [0; 1; 0; 0];

% 
Shaker_K = k_stinger*[-1 zeros(1, MDL.Ndofs-2) 1];

%% 
fsamp = 2^14;
T0 = 0;
T1 = 1;
dt = 1/fsamp;
Wfrc = 140;
Vamp = 2.8;

E = HARMONICSTIFFNESS(MDL.M, MDL.C, MDL.K, 2*pi*Wfrc, 1);
Uh = E\[Fshape*Vamp; Fshape*0];

U0 = Uh(1:MDL.Ndofs);
Ud0 = 2*pi*Wfrc*Uh(MDL.Ndofs+(1:MDL.Ndofs));

FEX = @(t) Fshape*Vamp*cos(2*pi*Wfrc*t);
% rng(1);
% FEX = Fshape.*rand(1, length(T0:dt:T1));

opts = struct('Display', 'waitbar');

tic 
[T, U, Ud, Udd, MDLh] = MDL.HHTAMARCH(T0, T1, dt, U0, Ud0, ...
                    FEX, opts);
toc

%%
figure(10)
clf()
% hold on
% plot(atan2(sin(2*pi*Wfrc*T), cos(2*pi*Wfrc*T)), U(1,:), '.')
plot(T, U(1,:), '.-')
ylabel('Displacement (m/s^2)')

figure(20)
clf()
% yyaxis right
% plot(atan2(sin(2*pi*Wfrc*T), cos(2*pi*Wfrc*T)), Shaker_K*U, '.')
plot(T, Shaker_K*U, '.-');
ylabel('Force (N)')
%%
% [freqs, Uf] = FFTFUN(T, U(1,:)');

% figure(1)
% clf()
% semilogy(freqs, abs(Uf))