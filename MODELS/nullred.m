clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM')

model = 'BRB';

load(sprintf('./%s/MATRICES.mat',model), 'M', 'K', 'Fv', 'R');
R = R';
Fv = Fv';
% Load Mesh
Nds = dlmread(sprintf('./%s/Nodes.dat', model));
Quad = dlmread(sprintf('./%s/Elements.dat', model));

MESH = MESH2D(Nds, 3, [], Quad, 2);

% Relative Coordinates [XT-XB; XB; eta]
Ngen = size(K, 1) - MESH.Nn*MESH.dpn*2;
Trel = [kron([1 1; 0 1], eye(MESH.Nn*MESH.dpn)), zeros(MESH.Nn*MESH.dpn*2,Ngen);
	zeros(Ngen, MESH.Nn*MESH.dpn*2), eye(Ngen)];

Mrel = Trel'*M*Trel;  Mrel = 0.5*(Mrel+Mrel');
Krel = Trel'*K*Trel;  Krel = 0.5*(Krel+Krel');
Rrel = R*Trel;
Fvrel = Trel'*Fv;

%% HCB Reduction
[Mhcb, Khcb, Th] = HCBREDUCE(Mrel, Krel, 1:MESH.Nn*MESH.dpn, Ngen);

Mhcb = 0.5*(Mhcb+Mhcb');
Khcb = 0.5*(Khcb+Khcb');

Rhcb = Rrel*Th;
Fvhcb = Th'*Fvrel;
Thcb = Trel*Th;

%% Fixed Interface RBMs
[V, D] = eigs(Khcb(MESH.Nn*MESH.dpn+1:end, MESH.Nn*MESH.dpn+1:end),
	      Mhcb(MESH.Nn*MESH.dpn+1:end, MESH.Nn*MESH.dpn+1:end), 10, 'SM');
[~, si] = sort(diag(D));
V = V(:, si);
L = null([zeros(MESH.Nn*MESH.dpn, 6); V(:, 1:6)]'*Mhcb);

K = L'*Khcb*L;
M = L'*Mhcb*L;
R = Rhcb*L;
Fv = L'*Fvhcb;
TFM = Thcb*L;  % Recover original system

disp('DONE!')

% save(sprintf('./%s/MATRICES_NR.mat',model), 'K', 'M', 'R', 'Fv', 'TFM', 'MESH', '-v7');
save(sprintf('./%s/MATRICES_NR.mat',model), 'K', 'M', 'R', 'L',  'Fv', 'TFM', '-v7');
