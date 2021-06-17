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
load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'K', 'M', 'R', 'L', 'Fv', 'TFM');

%% Static Natural Frequencies : Experimental
Wsex = zeros(3,1);  Zsex = zeros(3,1);
for i=1:3
    exp(1) = load(sprintf('./MATFILES/Mode%d_Low.mat',i), 'AMP_avg', 'FRE_avg', 'DAM_avg');
    exp(2) = load(sprintf('./MATFILES/Mode%d_Med.mat',i), 'AMP_avg', 'FRE_avg', 'DAM_avg');
    exp(3) = load(sprintf('./MATFILES/Mode%d_High.mat',i), 'AMP_avg', 'FRE_avg', 'DAM_avg');

    fav = robustfit([exp(1).AMP_avg; exp(2).AMP_avg; exp(3).AMP_avg], ...
                    [exp(1).FRE_avg; exp(2).FRE_avg; exp(3).FRE_avg]);
    zav = robustfit([exp(1).AMP_avg; exp(2).AMP_avg; exp(3).AMP_avg], ...
                    [exp(1).DAM_avg; exp(2).DAM_avg; exp(3).DAM_avg]);
    Wsex(i) = fav(1)*2*pi;
    Zsex(i) = zav(1);
end

%%
% Experimental Data (Mode 1)
exp(1) = load('./MATFILES/Mode1_Low.mat', 'AMP_avg', 'FRE_avg', 'DAM_avg');
exp(2) = load('./MATFILES/Mode1_Med.mat', 'AMP_avg', 'FRE_avg', 'DAM_avg');
exp(3) = load('./MATFILES/Mode1_High.mat', 'AMP_avg', 'FRE_avg', 'DAM_avg');

% ps = [0.05 0.50 0.99];
ps = [0.25 0.05];
ecs = colormap(lines(length(ps)));
% ecs = {'red', 'yellow', 'green'};
fcs = kron(ps(:)/2, [1 1 1]);
falph = 0.2;
pdists = {makedist('exp'), makedist('normal'), makedist('normal'), makedist('normal'), makedist('normal'), makedist('normal')};
polfuns = {@(ns, xs) PLAGU(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs)};
quadfuns = {@(n) LAGWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n)};
Pars = {'\mu', '\lambda', 'P', '\theta_X', '\theta_Y', 'gap'};
mdis = [1 3 5];
Nsamps = 100000;

is = {[1 2 3], [1 3 4]};
Nq_pce = {10};
Npce = 6;
calc = false;
for i=1:length(is)
    Nq_pces = ones(1, 6);
    Nq_pces(is{i}) = Nq_pce{i};

    Ir = cell(size(is{i}));
    [Ir{:}] = ndgrid(0:Nq_pce{i}-1);
    Ir = cell2mat(cellfun(@(c) c(:)', Ir(:), 'UniformOutput', false))';
    nxis = Ir*Nq_pce{i}.^((1:length(is{i}))'-1);

    Irr = zeros(Nq_pce{i}^length(is{i}), 6);
    Irr(:, is{i}) = Ir;

    % Construct PCE Coefficients
    pref = sprintf('nlbb_%s', sprintf('%d', is{i}));
    fname = sprintf('./ALLPCE/%s_cofs.mat', pref);
    try
        load(fname, 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'IJs', 'Integs')
        if calc
            error('calculate')
        end
    catch me
        IJs = Irr(sum(Irr, 2)<=Npce, :);  % Polynomial Order

        Xis = zeros(Nq_pce{i}, 6);
        Wis = zeros(Nq_pce{i}, 6);
        for j=1:6
            [Xis(:, j), Wis(:, j)] = quadfuns{j}(Nq_pces(j));
        end

        Rcofs     = zeros(100, size(IJs,1));
        Wcofs     = zeros(100, size(IJs,1));
        Zcofs     = zeros(100, size(IJs,1));
        Wstatcofs = zeros(10, size(IJs,1));
        for n=1:length(nxis)
            nxi = nxis(n);
            res = load(sprintf('./ALLPCE/%s/%s_%d.mat', pref, pref, nxi), 'Qs', 'Phi', 'Lams', 'Zts', 'Wstat');

            Psi = 1;
            wi = 1;
            Integs = 1;
            for j=1:6
                [psi, Intg] = polfuns{j}(IJs(:,j), Xis(Irr(n, j)+1, j));
                Psi = Psi.*psi;
                Integs = Integs.*Intg;

                wi = Wis(Irr(n, j)+1, j)*wi;
            end

            Qs = res.Qs;
            Rx = real((res.Qs.*(R(3,:)*res.Phi)').*res.Lams);
            Rcofs = Rcofs + wi*(Rx.*Psi)./Integs';
            Wcofs = Wcofs + wi*(sqrt(res.Lams).*Psi)./Integs';
            Zcofs = Zcofs + wi*(res.Zts.*Psi)./Integs';
            Wstatcofs = Wstatcofs + wi*(res.Wstat.*Psi)./Integs';

            fprintf('%d\n', n)
        end
        save(fname, 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'IJs', 'Integs')
    end
    
    %% Sampling
    xs = zeros(Nsamps, length(is{i}));
    Psixs = ones(Nsamps, size(IJs, 1));
    for j=1:length(is{i})
        xs(:, j) = random(pdists{is{i}(j)}, Nsamps, 1);
        Psixs = Psixs.*polfuns{j}(IJs(:,j), xs(:, j));
    end

    WPI = prctile(Wcofs*Psixs', [ps 1-ps]*100, 2);
    ZPI = prctile(Zcofs*Psixs', [ps 1-ps]*100, 2);
    WstatPI = prctile(Wstatcofs*Psixs', [ps 1-ps]*100, 2);
        
    %% Variances
    Wvar = (Wcofs(:,2:end).^2)*Integs(2:end);
    Zvar = (Zcofs(:,2:end).^2)*Integs(2:end);
    Wstatvar = (Wstatcofs(:,2:end).^2)*Integs(2:end);
    
    % First Order Sobol Indices
    i1s = [];
    W_S1 = zeros(100, length(is{i}));
    Z_S1 = zeros(100, length(is{i}));
    Wstat_S1 = zeros(10, length(is{i}));
    for j=1:length(is{i})
        i1s = [i1s find(all(IJs(2:end, setdiff(1:6, is{i}(j)))==0, 2))+1];
        W_S1(:, j)     = ((Wcofs(:, i1s(:,j)).^2)*Integs(i1s(:,j)))./Wvar;
        Z_S1(:, j)     = ((Zcofs(:, i1s(:,j)).^2)*Integs(i1s(:,j)))./Zvar;
        Wstat_S1(:, j) = ((Wstatcofs(:, i1s(:,j)).^2)*Integs(i1s(:,j)))./Wstatvar;
    end
    
    % Second Order Sobol Indices
    i2s = [];
    ns = factorial(length(is{i}))/factorial(2);
    W_S2 = zeros(100, ns);
    Z_S2 = zeros(100, ns);
    Wstat_S2 = zeros(10, ns);
    k = 1;
    j1j2s = [];
    for j1=1:length(is{i})
        for j2=j1+1:length(is{i})
            i2s = [i2s find(all(IJs(2:end, setdiff(1:6, is{i}([j1 j2])))==0, 2))+1];
            ik = setdiff(i2s, i1s(:));
            W_S2(:, k)     = ((Wcofs(:, ik).^2)*Integs(ik))./Wvar;
            Z_S2(:, k)     = ((Zcofs(:, ik).^2)*Integs(ik))./Zvar;
            Wstat_S2(:, k) = ((Wstatcofs(:, ik).^2)*Integs(ik))./Wstatvar;
            
            j1j2s = [j1j2s; is{i}([j1 j2])];
            k = k+1;
        end
    end
    
    % Frequency Figures
    figure((i-1)*10+1)
    clf()
    colos = DISTINGUISHABLE_COLORS(size(W_S1,2), 'k');
%     semilogx(Rcofs(:,1)/9.81, Wcofs(:,1)/2/pi, 'k.-'); hold on
    ff = gobjects(length(ps),1);
    for j=1:length(ps)
        ff(j) = fill([Rcofs(:,1); Rcofs(end:-1:1,1)]/9.81, ...
            [WPI(:, j); WPI(end:-1:1, j+length(ps))]/2/pi, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph); hold on
        legend(ff(j), sprintf('$%d^{th}-%d^{th}$ Percentiles', ps(j)*100, (1-ps(j))*100));
    end
    set(gca, 'XScale', 'log')
    aa = gobjects(size(W_S1,2), 1);
    for j=1:length(exp)
        aa(1) = semilogx(exp(j).AMP_avg, exp(j).FRE_avg, 'k-', 'LineWidth', 2);
    end
    legend(aa(1), 'Experimental Data')
    xlabel('Response Amplitude (g)')
    ylabel('Natural Frequency (Hz)')
    yyaxis right
    for k=1:size(W_S1,2)
        aa(1+k) = semilogx(Rcofs(:,1)/9.81, W_S1(:,k), '.-', 'Color', colos(k,:), 'LineWidth', 2);
        legend(aa(1+k), sprintf('$S_{%s}$', Pars{is{i}(k)}))
    end
    ll=legend([ff;aa], 'Location', 'west');
    set(ll, 'Visible', 'off')
    ylabel('First Order Sobol Indices')
    ylim([-0.1 1.1])
    xlim([min(Rcofs(:,1)/9.81) 150])
    
    figure((i-1)*10+2)
    clf()
    colos = DISTINGUISHABLE_COLORS(size(W_S2,2), 'k');
%     semilogx(Rcofs(:,1)/9.81, Wcofs(:,1)/2/pi, 'k.-'); hold on
    for j=1:length(ps)
        ff(j) = fill([Rcofs(:,1); Rcofs(end:-1:1,1)]/9.81, ...
            [WPI(:, j); WPI(end:-1:1, j+length(ps))]/2/pi, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph); hold on
        legend(ff(j), sprintf('$%d^{th}-%d^{th}$ Percentiles', ps(j)*100, (1-ps(j))*100));
    end
    set(gca, 'XScale', 'log')
    aa = gobjects(size(W_S2,2)+1, 1);
    for j=1:length(exp)
        aa(1)=semilogx(exp(j).AMP_avg, exp(j).FRE_avg, 'k-', 'LineWidth', 2);
    end
    legend(aa(1), 'Experimental Data')
    xlabel('Response Amplitude (g)')
    ylabel('Natural Frequency (Hz)')
    yyaxis right
    for k=1:size(W_S2,2)
        aa(1+k) = semilogx(Rcofs(:,1)/9.81, W_S2(:,k), '.-', 'Color', colos(k,:), 'LineWidth', 2);
        legend(aa(1+k), sprintf('$S_{%s,%s}$', Pars{j1j2s(k,1)}, Pars{j1j2s(k,2)}))
    end
    legend([ff;aa], 'Location', 'west')
    ylabel('Second Order Sobol Indices')
    ylim([-0.1 1.1])
    xlim([min(Rcofs(:,1)/9.81) 150])
    
    % Damping Figures
    figure((i-1)*10+5)
    clf()
    colos = DISTINGUISHABLE_COLORS(size(Z_S1,2), 'k');
%     semilogx(Rcofs(:,1)/9.81, Zcofs(:,1)*100, 'k.-'); hold on
    for j=1:length(ps)
        fill([Rcofs(:,1); Rcofs(end:-1:1,1)]/9.81, ...
            [ZPI(:, j); ZPI(end:-1:1, j+length(ps))]*100, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph); hold on
    end
    set(gca, 'XScale', 'log')
    for j=1:length(exp)
        semilogx(exp(j).AMP_avg, exp(j).DAM_avg*100, 'k-', 'LineWidth', 2);
    end
    xlabel('Response Amplitude (g)')
    ylabel('Natural Frequency (Hz)')
    yyaxis right
    for k=1:size(Z_S1,2)
        semilogx(Rcofs(:,1)/9.81, Z_S1(:,k), '.-', 'Color', colos(k,:), 'LineWidth', 2)
    end
    ylabel('First Order Sobol Indices')
    ylim([-0.1 1.1])
    xlim([min(Rcofs(:,1)/9.81) 150])
    
    figure((i-1)*10+6)
    clf()
    colos = DISTINGUISHABLE_COLORS(size(Z_S2,2), 'k');
%     semilogx(Rcofs(:,1)/9.81, Zcofs(:,1)*100, 'k.-'); hold on
    for j=1:length(ps)
        fill([Rcofs(:,1); Rcofs(end:-1:1,1)]/9.81, ...
            [ZPI(:, j); ZPI(end:-1:1, j+length(ps))]*100, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph); hold on
    end
    set(gca, 'XScale', 'log')
    for j=1:length(exp)
        semilogx(exp(j).AMP_avg, exp(j).DAM_avg*100, 'k-', 'LineWidth', 2);
    end
    xlabel('Response Amplitude (g)')
    ylabel('Natural Frequency (Hz)')
    yyaxis right
    for k=1:size(W_S2,2)
        semilogx(Rcofs(:,1)/9.81, Z_S2(:,k), '.-', 'Color', colos(k,:), 'LineWidth', 2)
    end
    ylabel('Second Order Sobol Indices')
    ylim([-0.1 1.1])
    xlim([min(Rcofs(:,1)/9.81) 150])
    
    % Modal Frequencies
    figure((i-1)*10+4)
    clf()
    colos = DISTINGUISHABLE_COLORS(size([Wstat_S1 Wstat_S2],2));
    bb = bar([Wstat_S1(mdis, :) Wstat_S2(mdis, :)]);
    for k=1:size(Wstat_S1,2)
        bb(k).FaceColor = colos(k,:);
        legend(bb(k), sprintf('$S_{%s}$', Pars{is{i}(k)}));
    end
    for k=size(Wstat_S1,2)+(1:size(Wstat_S2,2))
        bb(k).FaceColor = colos(k,:);
        legend(bb(k), sprintf('$S_{%s, %s}$', Pars{j1j2s(k-size(Wstat_S1,2),1)}, Pars{j1j2s(k-size(Wstat_S1,2),2)}));
    end
    set(gca, 'yscale', 'log')
    legend(bb, 'Location', 'best')
    xlabel('Bending Mode ID')
    ylabel('Sobol Indices')
    
    figure((i-1)*10+9)
    clf()
    for im=1:length(mdis)
        subplot(length(mdis), 1, im)
        histogram(rmoutliers(Psixs*Wstatcofs(mdis(im),:)'/2/pi), 'Normalization', 'pdf'); hold on
        ylabel('PDF')
        title(sprintf('Mode %d', im))
        ylim([0 2.5])
        plot(Wsex(im)/2/pi*[1 1], ylim, 'k-', 'LineWidth', 2)
    end
    xlabel('Natural Frequency (Hz)')
end