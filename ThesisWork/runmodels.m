clc
clear all

%% Mean Models
% QMIN = [-5.5 -6 -8.75];  QMAX = [-2.625 -3.15 -3.77];
QMIN = [-5.5 -6 -8.75];  QMAX = [-2.625 -3.2 -3.77];

% mdi = 1;  % Mode of Interest
% AMIN = -0.5;  AMAX = 2.5;
% rqnm_exprsurf_pinning;

% mdi = 2;  % Mode of Interest
% AMIN = 0;  AMAX = 3;
% rqnm_exprsurf_pinning;

mdi = 3;  % Mode of Interest
AMIN = -2;  AMAX = 3;
rqnm_exprsurf_pinning;

%%
% rqnm_exprsurf_gappce;
% 
% rqnm_exprsurf_mscpce;
% 
% rqnm_exprsurf_mupce;
% 
% rqnm_exprsurf_prespce;
% 
% rqnm_exprsurf_rotpce;
