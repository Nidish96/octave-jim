clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/QUADRATURE/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

model = 'BRB_Thesis';
load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'R');

%% Calculate PCE
Pars = {'\mu', '\lambda', 'P', '\theta_X', '\theta_Y', 'gap', 'rad'};
pdists = [makedist('exp') repmat(makedist('normal'), 1, 6)];
polfuns = [{@(ns, xs) PLAGU(ns, xs)}, repmat({@(ns, xs) PHERM(ns, xs)}, 1, 6)];
quadfuns = {@(n) LAGWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n)};
mdis = [1 3 5];
Nsamps = 100000;

apref = 'nlbb';
%is = {[1], [2], [3], [4 5], [6], [4 6], [7], [1 2 7], [1 2 3 4 6 7]}; % [1 3 4]
%Nq_pce = {10, 10, 10, 10, 10, 10, 10, 10, 5};
%Npce = 9;

is = {[1 2 3 4 6 7]}; % [1 3 4]
Nq_pce = {5};
Npce = 4;

for i=length(is)
    Nq_pces = ones(1, 7);
    Nq_pces(is{i}) = Nq_pce{i};

    Ir = cell(size(is{i}));
    [Ir{:}] = ndgrid(0:Nq_pce{i}-1);
    Ir = cell2mat(cellfun(@(c) c(:)', Ir(:), 'UniformOutput', false))';
    nxis = Ir*Nq_pce{i}.^((1:length(is{i}))'-1);

    Irr = zeros(Nq_pce{i}^length(is{i}), 7);
    Irr(:, is{i}) = Ir;

    % Construct PCE Coefficients
    pref = sprintf('nlbb_%s', sprintf('%d', is{i}));
    dpref = sprintf('%s_%s', apref, sprintf('%d', is{i}));
    fname = sprintf('./ALLPCE/%s_cofs.mat', dpref);
   
    IJs = Irr(sum(Irr, 2)<=Npce, :);  % Polynomial Order

    Xis = zeros(Nq_pce{i}, 7);
    Wis = zeros(Nq_pce{i}, 7);
    for j=1:7
        [Xis(:, j), Wis(:, j)] = quadfuns{j}(Nq_pces(j));
    end

    Rcofs     = zeros(100, size(IJs,1));
    Wcofs     = zeros(100, size(IJs,1));
    Zcofs     = zeros(100, size(IJs,1));
    Wstatcofs = zeros(10, size(IJs,1));
    
    QPs = struct('Wstats', zeros(10, length(nxis)), ...
        'Rs', zeros(100, length(nxis)), 'Ws', zeros(100, length(nxis)), ...
        'Zs', zeros(100, length(nxis)), ...
        'xxis', zeros(length(nxis), 7), 'wxis', zeros(length(nxis), 7));
    for n=1:length(nxis)
        nxi = nxis(n);
        load(sprintf('./ALLPCE/%s/%s_%d_m1.mat', dpref, pref, nxi), 'Qs', 'Phi', 'Lams', 'Zts', 'Wstat');
    
        Psi = ones(1, size(IJs,1));
        wi = 1;
        Integs = 1;
        for j=1:7
            [psi, Intg] = polfuns{j}(IJs(:,j), Xis(Irr(n, j)+1, j));
            Psi = Psi.*psi;
            Integs = Integs.*Intg;

            wi = Wis(Irr(n, j)+1, j)*wi;
            
            QPs.xxis(n, j) = Xis(Irr(n, j)+1, j);
            QPs.wxis(n, j) = Wis(Irr(n, j)+1, j);
        end

        Qs = Qs;
%         Rx = real((Qs.*(R(3,:)*Phi)').*Lams);  % Response Acceleration Amplitude (m)
        Rx = real((Qs.*(R(3,:)*Phi)'));  % Response Amplitude (m)
        Rcofs = Rcofs + wi*(Rx.*Psi)./Integs';
        Wcofs = Wcofs + wi*(sqrt(Lams).*Psi)./Integs';
        Zcofs = Zcofs + wi*(Zts.*Psi)./Integs';
        Wstatcofs = Wstatcofs + wi*(Wstat.*Psi)./Integs';
        
        QPs.Rs(:, n) = Rx;
        QPs.Ws(:, n) = sqrt(Lams);
        QPs.Zs(:, n) = Zts;
        QPs.Wstats(:, n) = Wstat;

        fprintf('%d\n', n)
    end
    
    %% Regression-based Estimates
    xxis = zeros(length(nxis), 7);
    wxis = zeros(length(nxis), 7);
    for n=1:length(nxis)
        for j=1:7
            xxis(n, j) = Xis(Irr(n, j)+1, j);
            wxis(n, j) = Wis(Irr(n, j)+1, j);
        end
    end
    wxis = prod(wxis, 2);
    Psi = ones(length(nxis), size(IJs, 1));
    for j=1:7
        Psi = Psi.*polfuns{j}(IJs(:,j), xxis(:,j));
    end
    
    Wcofs = zeros(size(QPs.Wstats,1), size(IJs,1));
    for n=1:size(QPs.Wstats,1)
        Wcofs(n, :) = wxis'*(QPs.Wstats(n ,:)'.*Psi);
    end

    % Regression
    Rregcofs = (Psi\QPs.Rs')';
    Wregcofs = (Psi\QPs.Ws')';
    Zregcofs = (Psi\QPs.Zs')';
    Wstatregcofs = (Psi\QPs.Wstats')';

    %% Save
    if length(is{i})==1
        x = sym('x');
        assume(x, 'real')
        rx = Rcofs*transpose(polfuns{is{i}}(0:Npce, x));
        wx = Wcofs*transpose(polfuns{is{i}}(0:Npce, x));
        zx = Zcofs*transpose(polfuns{is{i}}(0:Npce, x));
        wsx = Wstatcofs*transpose(polfuns{is{i}}(0:Npce, x));
        rxps = zeros(length(Qs), Npce+1);
        wxps = zeros(length(Qs), Npce+1);
        zxps = zeros(length(Qs), Npce+1);
        wsxps = zeros(10, Npce+1);
        for iq=1:length(Qs)
            rxps(iq, :) = sym2poly(rx(iq));
            wxps(iq, :) = sym2poly(wx(iq));
            zxps(iq, :) = sym2poly(zx(iq));
        end
        for iw=1:10
            wsxps(iw, :) = sym2poly(wsx(iw));
        end
        save(fname, 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'IJs', 'Integs', 'rxps', 'wxps', 'zxps', 'wsxps')
    else
        save(fname, 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'IJs', 'Integs')
    end

    %% PCE Simulations
    Rsim = EVALPCE(QPs.xxis, Rcofs, IJs, polfuns);
    Wsim = EVALPCE(QPs.xxis, Wcofs, IJs, polfuns);
    Zsim = EVALPCE(QPs.xxis, Zcofs, IJs, polfuns);
    Wstatsim = EVALPCE(QPs.xxis, Wstatcofs, IJs, polfuns);
    
    Rsimr = EVALPCE(QPs.xxis, Rregcofs, IJs, polfuns);
    Wsimr = EVALPCE(QPs.xxis, Wregcofs, IJs, polfuns);
    Zsimr = EVALPCE(QPs.xxis, Zregcofs, IJs, polfuns);
    Wstatsimr = EVALPCE(QPs.xxis, Wstatregcofs, IJs, polfuns);
    
    Wregsim = EVALPCE(xxis, Wstatregcofs, IJs, polfuns);
end
