% clc
clear all

pref = 'meanmodelbb1';

%% Mode 1
% RQNM_EXPRSURF_PCEFUN(ones(1, 7), 0, ones(1, 7), pref, 1, [-7.5 -4]);

%% Radius PCE
Nq_pces = [ones(1, 6) 10];
Ixs = ones(1, 7);
pref = 'testbb';
RQNM_EXPRSURF_PCEFUN(Ixs, 0, Nq_pces, pref, 1, [-7.5 -4]);