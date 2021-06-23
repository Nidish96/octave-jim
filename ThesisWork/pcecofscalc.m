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
pdists = {makedist('exp'), makedist('normal'), makedist('normal'), makedist('normal'), makedist('normal'), makedist('normal'), makedist('normal')};
polfuns = {@(ns, xs) PLAGU(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs)};
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
    Wstats = zeros(10, length(nxis));
    clf()
    for n=1:length(nxis)
        nxi = nxis(n);
        load(sprintf('./ALLPCE/%s/%s_%d_m1.mat', dpref, pref, nxi), 'Qs', 'Phi', 'Lams', 'Zts', 'Wstat');
        semilogx(Qs, sqrt(Lams)/2/pi', '.-'); hold on
%         plot(Wstat(1,:)/2/pi, '.-'); hold on
    
        Wstats(:, n) = Wstat;

        Psi = 1;
        wi = 1;
        Integs = 1;
        for j=1:7
            [psi, Intg] = polfuns{j}(IJs(:,j), Xis(Irr(n, j)+1, j));
            Psi = Psi.*psi;
            Integs = Integs.*Intg;

            wi = Wis(Irr(n, j)+1, j)*wi;
        end

        Qs = Qs;
%         Rx = real((Qs.*(R(3,:)*Phi)').*Lams);  % Response Acceleration Amplitude (m)
        Rx = real((Qs.*(R(3,:)*Phi)'));  % Response Amplitude (m)
        Rcofs = Rcofs + wi*(Rx.*Psi)./Integs';
        Wcofs = Wcofs + wi*(sqrt(Lams).*Psi)./Integs';
        Zcofs = Zcofs + wi*(Zts.*Psi)./Integs';
        Wstatcofs = Wstatcofs + wi*(Wstat.*Psi)./Integs';

        fprintf('%d\n', n)
    end
    
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
end
