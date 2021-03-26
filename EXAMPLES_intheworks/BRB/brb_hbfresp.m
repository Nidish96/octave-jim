clc
clear all
addpath('../../ROUTINES/')
addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/SOLVERS/')
addpath('../../ROUTINES/FEM')

%% Load System Matrices
load('./MATS/MATRICES_NR.mat', 'K', 'M', 'R', 'L', 'Fv');
M = single(M);  % Mass
K = single(K);  % Stiffness
R = single(R);  % Recovery (for selected nodes)
L = single(L);  % Null-Space Transform Matrix Xphys = L*Xnullred;
Fv = single(Fv);  % Bolt Force Vector - Applies 1N force through all 3 bolts at the same time - Scale appropriately for use

%% Read Interface Mesh Information
Nds = single(dlmread('./MATS/Nodes.dat'));
Quad = single(dlmread('./MATS/Elements.dat'));
Ne = size(Quad,1);

Nq = 2;  % sqrt(Number of quadrature points per element)
MESH = MESH2D(Nds, 3, [], Quad, Nq);
[Qm, Tm] = MESH.ND2QP();  % Qm -> Ne*Nq^2, Nn; Tm -> Nn, Ne*Nq^2

%% Fixed Interface Modes (Can be constructed analytically since we're just using a HCB model)
Nint = MESH.Nn*MESH.dpn;  % Interface Degrees of Freedom
Ndofs = size(K,1);
Vfi = [zeros(Nint, Ndofs-Nint); eye(Ndofs-Nint)];  % Mode Shapes
Wfi = sqrt(diag(K(Nint+1:end, Nint+1:end)));  % Frequencies

Zetas = [0.002; 0.004];  % Modal Damping (assumed)
ab = [1./(2*Wfi(1:length(Zetas))) Wfi(1:length(Zetas))/2]\Zetas;

C = ab(1)*M + ab(2)*K;  % Proportional Damping Matrix

%% Setup MDOF Object
MDL = MDOFGEN(M, K, C, L);

ELDRYFRICT