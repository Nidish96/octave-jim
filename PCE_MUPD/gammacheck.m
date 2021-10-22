clc
clear all
addpath('../ROUTINES/QUADRATURE/')

%% Check
[xi, wi] = LAGWT(5);

dist = makedist('gamma', 1, 1);
% dist = makedist('exp');

pdf(dist, xi)