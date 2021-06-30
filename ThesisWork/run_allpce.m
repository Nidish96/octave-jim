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
% is = [1; 3; 4];
% pref = "nlbb_134";

%% "Mu" PCE
is = 1;

%% Generate Parameter Space
pref = sprintf('nlbb_%s', sprintf('%d', is));

Nq_pces = ones(1, 7);
Nq_pces(is) = Nq_pce;

Ir = cell(size(is));
[Ir{:}] = ndgrid(0:Nq_pce-1);
Ir = cell2mat(cellfun(@(c) c(:)', Ir(:), 'UniformOutput', false))';
nxis = Ir*Nq_pce.^((1:length(is))'-1);

Irr = zeros(Nq_pce^length(is), 7);
Irr(:, is) = Ir;

for n=1:Nq_pce
%     worker(n) = batch(@RQNM_EXPRSURF_PCEFUN, 0, {Irr(n,:)+1, nxis(n), Nq_pces, pref, 1, [-7.5, -3], 1, 'quad'});
    RQNM_EXPRSURF_PCEFUN(Irr(n,:)+1, nxis(n), Nq_pces, pref, 1, [-7.5, -3], 1, 'quad');
end