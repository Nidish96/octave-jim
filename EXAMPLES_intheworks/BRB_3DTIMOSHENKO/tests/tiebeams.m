clc
clear all
addpath('../../ROUTINES/FEM/BEAMS')

%% Set Geometric Properties
Lbb = 11.8*25.4e-3;  % Body length
Lin = 120e-3;  % Interface Length
wdt = 25.4e-3;

%% Set Material Properties
pars.E = 2e11;
pars.G = 2e11/(2*(1+0.3));
pars.rho = 7800;

pars.A = wdt*wdt;
pars.Iy = wdt^4/12;
pars.Iz = wdt^4/12;

pars.k1 = 5.0/6;
pars.k2 = 5.0/6;
pars.k3 = 5.0/6;

parsint = pars;  % Half beam portion
parsint.A = pars.A/2;
parsint.Iy = pars.Iy/8;
parsint.Iz = pars.Iz/2;

%% Set Finite Element Discretization
% Nein = 8;  % Number of elements in the interface
% Nebb = fix(Nein*2.5);  % Number of elements in the rest of the body
Nebb = 20;

% Coordinates of nodes (neutral axis)
BM = struct('X', linspace(0, Lbb, Nebb+1)', ...
	    'Y', zeros(Nebb+1,1), ...
	    'Z', zeros(Nebb+1,1), ...
            'WY', wdt*ones(Nebb+1,1), ...
            'WZ', wdt*ones(Nebb+1,1));

BM_bot = struct('X', linspace(0, Lbb, Nebb+1)', ...
		'Y', zeros(Nebb+1,1), ...
		'Z', -wdt/4*ones(Nebb+1,1), ...
		'WY', wdt*ones(Nebb+1,1), ...
		'WZ', wdt/2*ones(Nebb+1,1));
BM_top = struct('X', linspace(0, Lbb, Nebb+1)', ...
		'Y', zeros(Nebb+1,1), ...
		'Z', wdt/4*ones(Nebb+1,1), ...
		'WY', wdt*ones(Nebb+1,1), ...
		'WZ', wdt/2*ones(Nebb+1,1));

BM_con = struct('X', [BM_bot.X; BM_top.X], ...
		'Y', [BM_bot.Y; BM_top.Y], ...
		'Z', [BM_bot.Z; BM_top.Z], ...
		'WY', [BM_bot.WY; BM_top.WY], ...
		'WZ', [BM_bot.WZ; BM_top.WZ]);

%% Analytical Solution (Free-Free)
kLs = [0 4.73 7.8532 10.9956 14.1371 17.2787]';
WAs = kLs.^2/range(BM.X)^2*sqrt(pars.E*pars.Iy/(pars.rho*pars.A));

%% Construct Matrices: "Target system"
[M, K] = CONSTRUCT(BM, pars);

%% Construct Matrices: Bottom and Top halves
[Mb, Kb] = CONSTRUCT(BM_bot, parsint);
[Mt, Kt] = CONSTRUCT(BM_top, parsint);

Lcon = [kron(eye(Nebb+1), [eye(3), [0 wdt/4 0;-wdt/4 0 0;0 0 0]; zeros(3), eye(3)]);
	kron(eye(Nebb+1), [eye(3), [0 -wdt/4 0;wdt/4 0 0;0 0 0]; zeros(3), eye(3)]);];

ctob = [eye(3), [0 wdt/4 0;-wdt/4 0 0;0 0 0]; zeros(3), eye(3)];
ctot = [eye(3), [0 -wdt/4 0;wdt/4 0 0;0 0 0]; zeros(3), eye(3)];
Lcon = [kron(eye(Nebb+1), eye(6));
	kron(eye(Nebb+1), ctot*inv(ctob))];

Mcon = Lcon'*blkdiag(Mb, Mt)*Lcon;
Kcon = Lcon'*blkdiag(Kb, Kt)*Lcon;

Ltot = Lcon;

%% Modal Analysis
[V, Ws] = eig(K, M);
[Ws, si] = sort(sqrt(abs(diag(Ws))));
V = V(:,si);

[Vcon, Wscon] = eig(Kcon, Mcon);
[Wscon, si] = sort(sqrt(abs(diag(Wscon))));
Vcon = Ltot*Vcon(:,si);

mi = 7;
sc = 1e-2;
figure(1)
clf()

subplot(1,2,1)
DEPICTBEAM_TM3D(diff(BM.X), BM.WY, BM.WZ, ...
		[BM.X BM.Y BM.Z], V(:, mi)*sc, 'b', 0.1, 2);
axis equal
grid on
xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')
title('Target')

subplot(1,2,2)
DEPICTBEAM_TM3D(diff(BM_bot.X), BM_bot.WY, BM_bot.WZ, ...
		[BM_bot.X BM_bot.Y BM_bot.Z], Vcon(1:(Nebb+1)*6, mi)*sc, 'b', 0.1, 2);
DEPICTBEAM_TM3D(diff(BM_top.X), BM_top.WY, BM_top.WZ, ...
		[BM_top.X BM_top.Y BM_top.Z], Vcon((Nebb+1)*6+1:end, mi)*sc, 'r', 0.1, 2);
axis equal
grid on
xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')
title('SubComp')

disp('Analytical Solution')
disp(WAs);
disp('Numerical Solution')
disp([Ws(1:10) Wscon(1:10)]);
