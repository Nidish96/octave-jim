clc
clear all
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/CONTACTMODELS')

% Load Matrices
load('../MODELS/BRB/MATRICES.mat', 'M', 'K', 'Fv')

% Load Mesh
Nds = dlmread('../MODELS/BRB/Nodes.dat');
Quad = dlmread('../MODELS/BRB/Elements.dat');

MESH = MESH2D(Nds, 3, [], Quad, 2);

% Test Out Elastic Dry Friction on the Mesh
Uxyn = repmat([1 2 3], MESH.Nn, 1);
Zxy = zeros(2, MESH.Ne*MESH.Nq^2);

Pars = [1; 2; 3; 0.1];
pA = repmat(eye(4), MESH.Ne*MESH.Nq^2, 1);

% Pars = repmat([1;2;3;0.1], 1, MESH.Ne*MESH.Nq^2);

U_n = [reshape(Uxyn', MESH.Nn*MESH.dpn, 1); zeros(size(K,1)-MESH.Nn*MESH.dpn,1)];
pU = speye(length(U_n));
pU(:, 10:end) = [];

MESH = MESH.SETCFUN(@(u, z, ud, P) ELDRYFRICT(u, z, ud, P, 0));
[F, Zxy, dFdU, dFdUd, dFdP] = MESH.CONTACTEVAL(pU'*U_n, Zxy, pU'*U_n*0, Pars, pA, pU);

% The above form should be directly applicable in HHTA_NONLIN_HYST

%% TODO:
% 1. BUILD STANDARD RESIDUE FUNCTIONS OUTSIDE OF THE MESH STRUCTURE for RQNM
% 2. BUILD STANDARD RESIDUE FUNCTIONS OUTSIDE OF THE MESH STRUCTURE for HBM
