clc
clear all
% Order of Parameters : [mu, msc, prestress, rotx, roty, gap]

%% Total Number of Quadrature Points
Nq_pce = 10;

%% "All" PCE
% is = (1:6)';

%% (i, j, k, ...)^th parameter (N-D) 
% is = [4; 5; 6];
% is = [1 2 3 4 6];
is = [1; 3; 4];
pref = "nlbb_134";

%% "Mu" PCE
% is = 1;
% pref = "nlbb_mupce";

%% Generate Parameter Space
Nq_pces = ones(1, 6);
Nq_pces(is) = Nq_pce;

Ir = cell(size(is));
[Ir{:}] = ndgrid(0:Nq_pce-1);
Ir = cell2mat(cellfun(@(c) c(:)', Ir(:), 'UniformOutput', false))';
nxis = Ir*Nq_pce.^((1:length(is))'-1);

Irr = zeros(Nq_pce^length(is), 6);
Irr(:, is) = Ir;

n = 443;
% [443 454 479]
for n=[444 455 480]
    RQNM_EXPRSURF_PCEFUN(Irr(n,:)+1, nxis(n), Nq_pces, pref);
%     STATRES_EXPRSURF_PCEFUN(Irr(n,:)+1, nxis(n), Nq_pces, pref);
end