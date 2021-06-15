clc
clear all
% Order of Parameters : [mu, msc, prestress, rotx, roty, gap]

%% Total Number of Quadrature Points
Nq_pce = 10;

%% "All" PCE
is = (1:6)';

%% (i, j, k, ...)^th parameter (N-D) 
is = [4; 5; 6];

%% "Mu" PCE
is = 1;
pref = "nlbb_mupce";

%% Generate Parameter Space
Nq_pces = ones(1, 6);
Nq_pces(is) = Nq_pce;

Ir = cell(size(is));
[Ir{:}] = ndgrid(0:Nq_pce-1);
Ir = cell2mat(cellfun(@(c) c(:)', Ir, 'UniformOutput', false))';
nxis = Ir*Nq_pce.^(is-1);

Irr = zeros(Nq_pce^length(is), 6);
Irr(:, is) = Ir;

for n=1:length(nxis)
%     RQNM_EXPRSURF_PCEFUN(Irr(n,:)+1, nxis(n), Nq_pces, pref);
    STATRES_EXPRSURF_PCEFUN(Irr(n,:)+1, nxis(n), Nq_pces, pref);
end