clc
clear all
addpath('../ROUTINES/HARMONIC/')

%% Parameters
% Linear Part
pars.m = 1.0;
pars.c = 0.2;
pars.k = 1.0;

% Bouc-Wen Parameters
pars.A  = 2.0;
pars.ki = 3.0;
pars.a  = 0.1;
pars.bt = 3.0;
pars.gm = 1.5;
pars.n  = 3;

%% Matrix Form 
M = [pars.m 0;0 0];
C = [pars.c 0;-pars.A 1];
K = [pars.k+pars.a*pars.ki (1-pars.a)*pars.ki;0 0];
E = [0; 1];  % "Shape" of non-linearity
