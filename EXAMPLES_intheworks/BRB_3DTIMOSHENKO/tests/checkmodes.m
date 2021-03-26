clc
clear all
addpath('../../ROUTINES/FEM/BEAMS')

%% Set Geometric Properties
Lbb = 11.8*25.4e-3;  % Body length
Lin = 120e-3;  % Interface Length
wdt = 25.4e-3;
wdt = 1e-3;

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
Nebb = 100;

% Coordinates of nodes (neutral axis)
BM1 = struct('X', linspace(0, Lbb, Nebb+1)', ...
	     'Y', zeros(Nebb+1,1), ...
	     'Z', zeros(Nebb+1,1), ...
             'WY', wdt*ones(Nebb+1,1), ...
             'WZ', wdt*ones(Nebb+1,1));

%% Analytical Solution (Free-Free)
kLs = [0 4.73 7.8532 10.9956 14.1371 17.2787]';
WAs = kLs.^2/range(BM1.X)^2*sqrt(pars.E*pars.Iy/(pars.rho*pars.A));

%% Construct Matrices
[M1, K1] = CONSTRUCT(BM1, pars);

%% Modal Analysis
[V, Ws] = eig(K1, M1);
[Ws, si] = sort(sqrt(abs(diag(Ws))));
V = V(:,si);

mi = 7;
sc = 1e-1;
figure(1)
clf()
DEPICTBEAM_TM3D(diff(BM1.X), BM1.WY, BM1.WZ, ...
		[BM1.X BM1.Y BM1.Z], V(:, mi)*sc, 'b', 0.1, 2);
axis equal
grid on
xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')

disp('Analytical Solution')
disp(WAs);
disp('Numerical Solution')
disp(Ws(1:10));
