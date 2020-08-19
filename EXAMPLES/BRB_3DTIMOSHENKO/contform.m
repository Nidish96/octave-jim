clc
clear all
addpath('../../ROUTINES/FEM')
addpath('../../ROUTINES/FEM/BEAMS')

%% Set Geometric Properties
Lbb = 11.8*25.4e-3;  % Body length
Lin = 120e-3;  % Interface Length
wdt = 25.4e-3;

%% Set Material Properties
pars.E = 2e11;
pars.rho = 7800;
pars.G = pars.E/(2*(1+0.3));

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

%% Analytical Solution
kLs = [0 4.73 7.8532 10.9956 14.1371 17.2787]';
WAs = kLs.^2/range([Beam1.X; Beam2.X])^2*sqrt(pars.E*pars.Iy/(pars.rho*pars.A))/2/pi;

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

dz = BM1.Z(1)-IN1.Z(end);
dy = BM1.Y(1)-IN1.Y(end);
L1(bis, eis) = [eye(3), [0 dz -dy; -dz 0 0; dy 0 0];
		zeros(3), eye(3)];  % Represent body node as a function of interface node
L1(:, bis) = [];   % Only retain node on the interface

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

dz = BM2.Z(1)-IN2.Z(end);
dy = BM2.Y(1)-IN2.Y(end);
L2(bis, eis) = [eye(3), [0 -dz dy; dz 0 0; -dy 0 0];
		zeros(3), eye(3)];   % Represent body node as a function of interface node
L2(:, bis) = [];   % Retain node on the interface

Mbm2 = L2'*blkdiag(Mi2, M2)*L2;
Kbm2 = L2'*blkdiag(Ki2, K2)*L2;

%% SET UP QUADRATURE
No = 2;

Les = diff(IN1.X);
Wys = IN1.WY(1:end-1);
Zs  = 0-IN1.Z(1:end-1);
[Q1, T1] = TM3D_ND2QP(Les, Wys, Zs, No);

Les = diff(IN2.X);
Wys = IN2.WY(1:end-1);
Zs  = 0-IN2.Z(1:end-1);
[Q2, T2] = TM3D_ND2QP(Les, Wys, Zs, No);

% Assembly for relative coordinates at quadrature locations
int1nds = Nebb+(1:(Nein+1));
int2nds = (Nebb+Nein+1)+(1:(Nein+1));
remnds  = setdiff(1:(Nebb+Nein+1)*2, [int1nds, int2nds]);

i1s = reshape((int1nds-1)*6+(1:6)', [], 1);
i2s = reshape((int2nds-1)*6+(1:6)', [], 1);

%% IN1 is facing downwards (z-) and IN2 is facing upwards (z+) so
%% positive interference is detected by IN2_disp-IN1_disp and this is
%% the sense in which relative displacements will be interpreted.

Qrel = zeros(Nein*No^2*3, (Nebb+Nein+1)*2*6);
Trel = zeros((Nebb+Nein+1)*2*6, Nein*No^2*3);

Qrel(:, i1s) = -Q1;
Qrel(:, i2s) = Q2;

Trel(i1s, :) = -T1;
Trel(i2s, :) = T2;

%% Recovery Matrix
recovnds = [1 (Nebb+Nein+1)*2];  % First and last nodes

%% Full System Assembly
M = blkdiag(Mbm1, Mbm2);
K = blkdiag(Kbm1, Kbm2);

%% Contact Model (Hertz-Mindlin)
Prestress = 12e3;
Aint = sum(sum(T1(1:6:end, :)));

nu   = 0.3;
Pint = Prestress/Aint;
sint = 10e-6;
chi  = 2.0;
ktkn = chi*(1-nu)/(2-nu);
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn   = kt/ktkn;

tstiff = zeros(Nein*No^2*3,1);
tstiff(1:3:end) = kt;
tstiff(2:3:end) = kt;
tstiff(3:3:end) = kn;

Ks = Trel*diag(tstiff)*Qrel;

[V, W] = eig(K, M);
[W, si] = sort(sqrt(diag(W))/2/pi);
V = V(:,si);
V = V./sqrt(diag(V'*M*V))';
V1 = V(1:(Nebb+Nein+1)*6, :);
V2 = V((Nebb+Nein+1)*6+1:end, :);

[V0, W0] = eig(K+Ks, M);
[W0, si] = sort(sqrt(diag(W0))/2/pi);
V0 = V0(:,si);
V0 = V0./sqrt(diag(V0'*M*V0))';
V01 = V0(1:(Nebb+Nein+1)*6, :);
V02 = V0((Nebb+Nein+1)*6+1:end, :);

%% Plot Mode Shapes

mi = 13;
sc = 1e-2;

figure(1)
clf()
DEPICTBEAM_TM3D(diff(Beam1.X), Beam1.WY, Beam1.WZ, ...
		[Beam1.X, Beam1.Y, Beam1.Z], L1*V01(:, mi)*sc, 'b',
		0.1, 2);
DEPICTBEAM_TM3D(diff(Beam2.X), Beam2.WY, Beam2.WZ, ...
		[Beam2.X, Beam2.Y, Beam2.Z], L2*V02(:, mi)*sc, 'r',
		0.1, 2); 

title(sprintf('Mode %d: %.2f Hz', mi, W0(mi)))
axis equal
grid on
xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')
