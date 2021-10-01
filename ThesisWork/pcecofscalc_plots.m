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
quadfuns = {@(n) LAGWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n ), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n)};
mdis = [1 3 5];
Nsamps = 100000;

apref = 'nlbb';
% is = {[1], [2], [3], [4 5], [6], [4 6], [7], [1 2 7], [1 2 3 4 6 7]}; % [1 3 4]
% Nq_pce = {10, 10, 10, 10, 10, 10, 10, 10, 5};
% Npce = 4;

is = {[1], [2], [3], [4 5], [6], [7]}; % [1 3 4]
Nq_pce = [10, 10, 10, 10, 10, 10, 10, 10, 5];
Npce = 4;

for i=4
    Nq_pces = ones(1, 7);
    Nq_pces(is{i}) = Nq_pce(i);

    Ir = cell(size(is{i}));
    [Ir{:}] = ndgrid(0:Nq_pce(i)-1);
    Ir = cell2mat(cellfun(@(c) c(:)', Ir(:), 'UniformOutput', false))';
    nxis = Ir*Nq_pce(i).^((1:length(is{i}))'-1);

    Irr = zeros(Nq_pce(i)^length(is{i}), 7);
    Irr(:, is{i}) = Ir;

    % Construct PCE Coefficients
    pref = sprintf('nlbb_%s', sprintf('%d', is{i}));
    dpref = sprintf('%s_%s', apref, sprintf('%d', is{i}));
    fname = sprintf('./ALLPCE/%s_cofs.mat', dpref);
   
    IJs = Irr(sum(Irr, 2)<=Npce, :);  % Polynomial Order

    Xis = zeros(Nq_pce(i), 7);
    Wis = zeros(Nq_pce(i), 7);
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

    %% Basis Functions At Evaluation Points
    Psi = ones(length(nxis), size(IJs, 1));
    for j=1:7
        Psi = Psi.*polfuns{j}(IJs(:,j), QPs.xxis(:,j));
    end

    %% POD-PCE 
    rm = mean(QPs.Rs,2); wm = mean(QPs.Ws,2); zm = mean(QPs.Zs,2); 
    [Ur, Sr, Vr] = svd(QPs.Rs-rm, 'econ');
    [Uw, Sw, Vw] = svd(QPs.Ws-wm, 'econ');
    [Uz, Sz, Vz] = svd(QPs.Zs-zm, 'econ');
%     [U, S, V] = svd([(QPs.Rs-rm)/mean(rm); (QPs.Ws-wm)/mean(wm); (QPs.Zs-zm)/mean(zm)], 'econ');
    [U, S, V] = svd([QPs.Rs./rm; QPs.Ws./wm; QPs.Zs./zm], 'econ');

%     figure(1)
%     clf()
%     semilogy(diag(Sr)/Sr(1,1), '.-'); hold on
%     semilogy(diag(Sw)/Sw(1,1), '.-')
%     semilogy(diag(Sz)/Sz(1,1), '.-')
%     semilogy(diag(S)/S(1,1), '.-')
%     xlabel('Index')
%     ylabel('Relative Singular Value')

    thresh = 1e-3;
    k = find(diag(Sr)/Sr(1,1)<thresh & diag(Sw)/Sw(1,1)<thresh & diag(Sz)/Sz(1,1)<thresh, 1);
    if isempty(k)
        k = size(S,1);
    end
    POD.R = U(1:100, 1:k).*rm;
    POD.W = U(101:200, 1:k).*wm;
    POD.Z = U(201:300, 1:k).*zm;
    POD.as = S(1:k, 1:k)*V(:, 1:k)';
    POD.acofs = zeros(k, size(IJs,1));
%     for j=1:k
%         POD.acofs(j, :) = prod(QPs.wxis,2)'*(POD.as(j,:)'.*Psi);
%     end
    POD.acofs = (Psi\POD.as')';
    
    %% Regression-based Estimates
    % Regression
    Rregcofs = (Psi\QPs.Rs')';
    Wregcofs = (Psi\QPs.Ws')';
    Zregcofs = (Psi\QPs.Zs')';
    Wstatregcofs = (Psi\QPs.Wstats')';
    
%     if i==1
%         Rregcofs = (Psi(2:end,:)\QPs.Rs(:,2:end)')';
%         Wregcofs = (Psi(2:end,:)\QPs.Ws(:,2:end)')';
%         Zregcofs = (Psi(2:end,:)\QPs.Zs(:,2:end)')';
%         Wstatregcofs = (Psi(2:end,:)\QPs.Wstats(:,2:end)')';
%     end
    
    %% x points for Simulations
    Nx = 100;
    xis = zeros(Nx, length(is{i}));
    for j=1:length(is{i})
        if is{i}(j)==1
            xis(:, j) = linspace(0, 31, Nx);
        else
            xis(:, j) = linspace(-1, 1, Nx)*6;
        end
    end
    switch length(is{i})
        case 1
            xis = xis;
        case 2
            [xxis, etais] = meshgrid(xis(:, 1), xis(:, 2));
            xis = [xxis(:) etais(:)];
        otherwise
            error('only 2D maximum')
    end
    xs = zeros(size(xis,1), 7);
    xs(:, is{i}) = xis;    
    
    %% PCE Simulations
    % Projection-Based PCE
%     Rsim = EVALPCE(QPs.xxis, Rcofs, IJs, polfuns);
%     Wsim = EVALPCE(QPs.xxis, Wcofs, IJs, polfuns);
%     Zsim = EVALPCE(QPs.xxis, Zcofs, IJs, polfuns);
    Wstatsim = EVALPCE(QPs.xxis, Wstatcofs, IJs, polfuns);
    Wstatsim_x = EVALPCE(xs, Wstatcofs, IJs, polfuns);
    
    % Regression-based PCE
%     Rsimr = EVALPCE(QPs.xxis, Rregcofs, IJs, polfuns);
%     Wsimr = EVALPCE(QPs.xxis, Wregcofs, IJs, polfuns);
%     Zsimr = EVALPCE(QPs.xxis, Zregcofs, IJs, polfuns);
    Wstatsimr = EVALPCE(QPs.xxis, Wstatregcofs, IJs, polfuns);
    Wstatsimr_x = EVALPCE(xs, Wstatregcofs, IJs, polfuns);
    
%     % POD-based PCE
%     asim = EVALPCE(QPs.xxis, POD.acofs, IJs, polfuns);
%     Rsimp = POD.R*asim;
%     Wsimp = POD.W*asim;
%     Zsimp = POD.Z*asim;

    %% Plot Fits
    switch length(is{i})
        case 1
            figure(i)
            clf()
            plot(QPs.xxis(:,is{i}), QPs.Wstats(1,:)/2/pi, 'ko', 'MarkerFaceColor', 'k'); hold on
%             plot(QPs.xxis(2:end,is{i}), QPs.Wstats(1,2:end)/2/pi, 'ko', 'MarkerFaceColor', 'k'); hold on
            plot(xis, Wstatsim_x(1,:)/2/pi, 'b-', 'LineWidth', 2);
            plot(xis, Wstatsimr_x(1,:)/2/pi, 'g-', 'LineWidth', 2);
            xlabel(sprintf('Germ %d', is{i}))
            ylabel('Mode 1 Natural Frequency (Hz)')
            grid on
            ylim((mean(QPs.Wstats(1,:))+std(QPs.Wstats(1,:))*3*[-1 1])/2/pi)
            yyaxis right;
            plot(xis, pdf(pdists(is{i}), xis), 'r-')
            ylabel('Germ PDF')
            ax = gca;
            ax.YAxis(2).Color = 'r';
        case 2
            figure(i)
            clf()
            plot3(QPs.xxis(:,is{i}(1)), QPs.xxis(:,is{i}(2)), QPs.Wstats(1,:)/2/pi, 'ko', 'MarkerFaceColor', 'k'); hold on
%             surf(reshape(xis(:,1), Nx, Nx), reshape(xis(:,2), Nx, Nx), reshape(Wstatsimr_x(1,:)/2/pi, Nx, Nx), 'EdgeColor', 'none');
%             surf(reshape(xis(:,1), Nx, Nx), reshape(xis(:,2), Nx, Nx), reshape(Wstatsim_x(1,:)/2/pi, Nx, Nx), 'EdgeColor', 'none');

            surf(reshape(xis(:,1), Nx, Nx), reshape(xis(:,2), Nx, Nx), reshape(Wstatsimr_x(1,:)/2/pi, Nx, Nx), reshape(pdf(pdists(is{i}(1)), xis(:,1)).*pdf(pdists(is{i}(2)), xis(:,2)),Nx,Nx), 'EdgeColor', 'none');
            colormap(jet)
            yy=colorbar;
            ylabel(yy, 'Germ PDF', 'interpreter', 'latex')
            xlabel(sprintf('Germ %d', is{i}(1)))
            ylabel(sprintf('Germ %d', is{i}(2)))
            zlabel('Mode 1 Natural Frequency (Hz)')
            grid on
        otherwise
            error('only 2D maximum')
    end
    switch(i)
        case 1
            legend('Simulations', 'PCE Projection Fit', 'Regression Fit', 'Location', 'best')
        case 4
            legend('Simulations', 'PCE Fit', 'Location', 'best')
    end
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/PCESINGFIT_%d.eps', i), '-depsc')
end
