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
Nein = 8;  % Number of elements in the interface
Nebb = fix(Nein*2.5);  % Number of elements in the rest of the body

% Coordinates of nodes (neutral axis)
BM1 = struct('X', linspace(0, Lbb, Nebb+1)', ...
	     'Y', zeros(Nebb+1,1), ...
	     'Z', zeros(Nebb+1,1), ...
             'WY', wdt*ones(Nebb+1,1), ...
             'WZ', wdt*ones(Nebb+1,1));
IN1 = struct('X', linspace(Lbb, Lbb+Lin, Nein+1)', ...
	     'Y', zeros(Nein+1,1), ...
	     'Z', ones(Nein+1,1)*(-wdt/4), ...
             'WY', wdt*ones(Nein+1,1), ...
             'WZ', wdt/2*ones(Nein+1,1));
Beam1 = struct('X', [BM1.X; IN1.X], ...
	       'Y', [BM1.Y; IN1.Y], ...
	       'Z', [BM1.Z; IN1.Z], ...
               'WY', [BM1.WY; IN1.WY], ...
               'WZ', [BM1.WZ(1:end-1); 0; IN1.WZ]);

IN2 = struct('X', linspace(Lbb, Lbb+Lin, Nein+1)', ...
	     'Y', zeros(Nein+1,1), ...
	     'Z', ones(Nein+1,1)*wdt/4, ...
             'WY', wdt*ones(Nein+1,1), ...
             'WZ', wdt/2*ones(Nein+1,1));
BM2 = struct('X', linspace(Lbb+Lin, Lbb+Lin+Lbb, Nebb+1)', ...
	     'Y', zeros(Nebb+1,1), ...
	     'Z', zeros(Nebb+1,1), ...
             'WY', wdt*ones(Nebb+1,1), ...
             'WZ', wdt*ones(Nebb+1,1));
Beam2 = struct('X', [IN2.X; BM2.X], ...
	       'Y', [IN2.Y; BM2.Y], ...
	       'Z', [IN2.Z; BM2.Z], ...
               'WY', [IN2.WY; BM2.WY], ...
	       'WZ', [IN2.WZ(1:end-1); 0; BM2.WZ]);

%% Construct Matrices

%% BEAM 1
% Body 1
[M1, K1] = CONSTRUCT(BM1, pars);
% Interface 1
[Mi1, Ki1] = CONSTRUCT(IN1, parsint);
% Tie the last node of BM1 to 1st node of IN1 by placing IN1 on BM1's section
L1 = eye((Nebb+1+Nein+1)*6);
bnode = Nebb+1;  % Last node of body
bis   = (bnode-1)*6+(1:6);
inode = Nebb+1+1;  % First node in interface
eis   = (inode-1)*6+(1:6);

% dz = IN1.Z(1)-BM1.Z(end);
% dy = IN1.Y(1)-BM1.Y(end);
% L1(eis, bis) = [eye(3), [0 -dz dy; dz 0 0; -dy 0 0];
% 		zeros(3), eye(3)];
% L1(:, eis) = [];

dz = BM1.Z(1)-IN1.Z(end);
dy = BM1.Y(1)-IN1.Y(end);
L1(bis, eis) = [eye(3), [0 -dz dy; dz 0 0; -dy 0 0];
		zeros(3), eye(3)];
L1(:, bis) = [];

Mbm1 = L1'*blkdiag(M1, Mi1)*L1;
Kbm1 = L1'*blkdiag(K1, Ki1)*L1;

%% BEAM 2
% Interface 2
[Mi2, Ki2] = CONSTRUCT(IN2, parsint);
% Body 2
[M2, K2] = CONSTRUCT(BM2, pars);
% Tie the last node of BM1 to 1st node of IN1 by placing IN1 on BM1's section
L2 = eye((Nein+1+Nebb+1)*6);
inode = Nein+1;  % Last node of interface
eis   = (inode-1)*6+(1:6);
bnode = Nein+1+1;  % First node of body
bis   = (bnode-1)*6+(1:6);

% dz = IN2.Z(1)-BM2.Z(end);
% dy = IN2.Y(1)-BM2.Y(end);
% L2(eis, bis) = [eye(3), [0 -dz dy; dz 0 0; -dy 0 0];
% 		zeros(3), eye(3)];
% L2(:, eis) = [];

dz = BM2.Z(1)-IN2.Z(end);
dy = BM2.Y(1)-IN2.Y(end);
L2(bis, eis) = [eye(3), [0 -dz dy; dz 0 0; -dy 0 0];
		zeros(3), eye(3)];
L2(:, bis) = [];

Mbm2 = L2'*blkdiag(Mi2, M2)*L2;
Kbm2 = L2'*blkdiag(Ki2, K2)*L2;


%% Check Mode shapes
[V1, Ws1] = eig(Kbm1, Mbm1);
[Ws1, si] = sort(sqrt(abs(diag(Ws1)))/2/pi);
V1 = V1(:, si);

[V2, Ws2] = eig(Kbm2, Mbm2);
[Ws2, si] = sort(sqrt(abs(diag(Ws2)))/2/pi);
V2 = V2(:, si);

V1 = L1*V1;
V2 = L2*V2;

[Ws1(1:10) Ws2(1:10)]

mi = 7;
sc = 2.5e-1;
sc = 5e-2;
figure(1)
clf()
DEPICTBEAM_TM3D(diff(Beam1.X), Beam1.WY, Beam1.WZ, ...
		[Beam1.X Beam1.Y Beam1.Z], V1(:, mi)*sc, 'b', 0.1, 2);
% DEPICTBEAM_TM3D(diff(Beam2.X), Beam2.WY, Beam2.WZ, ...
% 		[Beam2.X Beam2.Y Beam2.Z], V2(:, mi)*sc, 'r', 0.1, 2);

grid on
axis equal
title(sprintf('Mode %d: %.2f Hz', mi, Ws1(mi)))
xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')

print('./md15.eps', '-depsc')
