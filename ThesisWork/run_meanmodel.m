% clc
clear all

% pref = 'meanmodelbb1';

%% Mode 1
pref = 'FLATMEANMODELS/mmbb';
RQNM_EXPRSURF_PCEFUN(ones(1, 7), 0, ones(1, 7), pref, 1, [-7.5 -3], 0);

%% Mode 2
pref = 'FLATMEANMODELS/mmbb';
RQNM_EXPRSURF_PCEFUN(ones(1, 7), 0, ones(1, 7), pref, 2, [-8.5 -4.5], 0);

%% Mode 3
pref = 'FLATMEANMODELS/mmbb';
RQNM_EXPRSURF_PCEFUN(ones(1, 7), 0, ones(1, 7), pref, 3, [-9.5 -4.5], 0);
