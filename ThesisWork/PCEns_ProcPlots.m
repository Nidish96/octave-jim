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
    
    for j=1:3
        ni = find(all(isfinite([exp(j).AMP_avg exp(j).FRE_avg exp(j).DAM_avg]), 2));
        exp(j).AMP_avg = exp(j).AMP_avg(ni);
        exp(j).FRE_avg = exp(j).FRE_avg(ni);
        exp(j).DAM_avg = exp(j).DAM_avg(ni);
    end

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

for j=1:3
    ni = find(all(isfinite([exp(j).AMP_avg exp(j).FRE_avg exp(j).DAM_avg]), 2));
    exp(j).AMP_avg = exp(j).AMP_avg(ni);
    exp(j).FRE_avg = exp(j).FRE_avg(ni);
    exp(j).DAM_avg = exp(j).DAM_avg(ni);
%     exp(i).AMP_avg = exp(i).AMP_avg/9.81;
    exp(j).AMP_avg = exp(j).AMP_avg./(2*pi*exp(j).FRE_avg).^2;
end

% ps = [0.05 0.50 0.99];
ps = [0.05 0.25 0.45];
% ecs = colormap(lines(length(ps)));
ecs = zeros(length(ps), 3);
lst = {'-', '--', '-.'};
% ecs = {'red', 'yellow', 'green'};
fcs = kron((1-ps(:)), [1 1 1]);
% ecs = fcs;
falph = 1;
pdists = {makedist('exp'), makedist('normal'), makedist('normal'), makedist('normal'), makedist('normal'), makedist('normal'), makedist('normal')};
polfuns = {@(ns, xs) PLAGU(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs), @(ns, xs) PHERM(ns, xs)};
quadfuns = {@(n) LAGWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n), @(n) GPHWT(n)};
Pars = {'\mu', '\lambda', 'P', '\theta_X', '\theta_Y', 'gap', 'rad'};
mdis = [1 3 5];
Nsamps = 1000000;

% is = {[1 2 7]};
% Nq_pce = {10};
is = {[1 2 3 4 6 7]};
Nq_pce = {5};
zeta0 = 1.3841e-4;
for i=1:length(is)
    Nq_pces = ones(1, length(Pars));
    Nq_pces(is{i}) = Nq_pce{i};

    Ir = cell(size(is{i}));
    [Ir{:}] = ndgrid(0:Nq_pce{i}-1);
    Ir = cell2mat(cellfun(@(c) c(:)', Ir(:), 'UniformOutput', false))';
    nxis = Ir*Nq_pce{i}.^((1:length(is{i}))'-1);

    Irr = zeros(Nq_pce{i}^length(is{i}), length(Pars));
    Irr(:, is{i}) = Ir;

    % Construct PCE Coefficients
    pref = sprintf('nlbb_%s', sprintf('%d', is{i}));
    fname = sprintf('./ALLPCE/%s_cofs.mat', pref);

    load(fname, 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'IJs', 'Integs')    
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
        
    %% Variances & Sobol Indices
    Wvar = (Wcofs(:,2:end).^2)*Integs(2:end);
    Zvar = (Zcofs(:,2:end).^2)*Integs(2:end);
    Wstatvar = (Wstatcofs(:,2:end).^2)*Integs(2:end);
    
    % First Order Sobol Indices
    i1s = [];
    W_S1 = zeros(100, length(is{i}));
    Z_S1 = zeros(100, length(is{i}));
    Wstat_S1 = zeros(10, length(is{i}));
    for j=1:length(is{i})
        i1s = [i1s find(all(IJs(2:end, setdiff(1:length(Pars), is{i}(j)))==0, 2))+1];
        W_S1(:, j)     = ((Wcofs(:, i1s(:,j)).^2)*Integs(i1s(:,j)))./Wvar;
        Z_S1(:, j)     = ((Zcofs(:, i1s(:,j)).^2)*Integs(i1s(:,j)))./Zvar;
        Wstat_S1(:, j) = ((Wstatcofs(:, i1s(:,j)).^2)*Integs(i1s(:,j)))./Wstatvar;
    end
    
    % Second Order Sobol Indices
    i2s = [];
    ns = factorial(length(is{i}))/(factorial(length(is{i})-2)*factorial(2));
    W_S2 = zeros(100, ns);
    Z_S2 = zeros(100, ns);
    Wstat_S2 = zeros(10, ns);
    k = 1;
    j1j2s = [];
    for j1=1:length(is{i})
        for j2=j1+1:length(is{i})
            i2s = [i2s find(all(IJs(2:end, setdiff(1:length(Pars), is{i}([j1 j2])))==0, 2))+1];
            ik = setdiff(i2s, i1s(:));
            W_S2(:, k)     = ((Wcofs(:, ik).^2)*Integs(ik))./Wvar;
            Z_S2(:, k)     = ((Zcofs(:, ik).^2)*Integs(ik))./Zvar;
            Wstat_S2(:, k) = ((Wstatcofs(:, ik).^2)*Integs(ik))./Wstatvar;
            
            j1j2s = [j1j2s; is{i}([j1 j2])];
            k = k+1;
        end
    end
    
    % Total Effect Sobol Indices
    iT1s = [];
    W_ST1 = zeros(100, length(is{i}));
    Z_ST1 = zeros(100, length(is{i}));
    Wstat_ST1 = zeros(10, length(is{i}));
    for j=1:length(is{i})
        iT1s = [iT1s find(IJs(2:end, is{i}(j))~=0)+1];
        W_ST1(:, j)     = ((Wcofs(:, iT1s(:,j)).^2)*Integs(iT1s(:,j)))./Wvar;
        Z_ST1(:, j)     = ((Zcofs(:, iT1s(:,j)).^2)*Integs(iT1s(:,j)))./Zvar;
        Wstat_ST1(:, j) = ((Wstatcofs(:, iT1s(:,j)).^2)*Integs(iT1s(:,j)))./Wstatvar;
    end
    
    W_S1 = sqrt(W_S1);
    W_S2 = sqrt(W_S2);
    W_ST1 = sqrt(W_ST1);
    Wstat_S1 = sqrt(Wstat_S1);
    Wstat_S2 = sqrt(Wstat_S2);
    Wstat_ST1 = sqrt(Wstat_ST1);
    %% S1 Plots    
    % Frequency - S1
    figure((i-1)*10+1)
    clf()
    set(gcf, 'Color', 'white')
    colos = DISTINGUISHABLE_COLORS(size(W_S1,2));
    ff = gobjects(length(ps),1);
    for j=1:length(ps)
        ff(j) = fill([Rcofs(:,1); Rcofs(end:-1:1,1)], ...
            [WPI(:, j); WPI(end:-1:1, j+length(ps))]/2/pi, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
        legend(ff(j), sprintf('$%d^{th}-%d^{th}$ Percentiles', fix(ps(j)*100), fix((1-ps(j))*100)));
    end
    set(gca, 'XScale', 'log')
    aa = gobjects(size(W_S1,2), 1);
    for j=1:length(exp)
        aa(1) = semilogx(exp(j).AMP_avg, exp(j).FRE_avg, 'k-', 'LineWidth', 2);
    end
    legend(aa(1), 'Experimental Data')
    xlabel('Response Amplitude (m)')
    ylabel('Natural Frequency (Hz)')
    
    ll=legend([ff;aa(1)], 'Location', 'southwest'); set(ll, 'Visible', 'off')
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    
    export_fig('./FIGS/ALLPCE_BBW_S1s_0.png', '-dpng', '-r300')
    
    yyaxis right
    for k=1:size(W_S1,2)
        aa(1+k) = semilogx(Rcofs(:,1), W_S1(:,k), '-', 'Color', colos(k,:), 'LineWidth', 2);
        legend(aa(1+k), sprintf('$S_{%s}$', Pars{is{i}(k)}))
    end
    ll=legend([ff;aa], 'Location', 'southwest');
    set(ll, 'Visible', 'off');
    ylabel('First Order Sobol Indices (sq. root)')
    ylim([-0.1 1.1])
%     xlim(minmax(Rcofs(:,1)'))
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_BBW_S1s.png', '-dpng', '-r300')
    
    %% Damping - S1
    figure((i-1)*10+5)
    clf()
    set(gcf, 'Color', 'white')
    for j=1:length(ps)
        fill([Rcofs(:,1); Rcofs(end:-1:1,1)], ...
            [ZPI(:, j); ZPI(end:-1:1, j+length(ps))]*100+zeta0*100, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
    end
    set(gca, 'XScale', 'log')
    for j=1:length(exp)
        semilogx(exp(j).AMP_avg, exp(j).DAM_avg*100, 'k-', 'LineWidth', 2);
    end
    xlabel('Response Amplitude (m)')
    ylabel('Damping Factor (\%)')
    yl = ylim;
    ylim([0 yl(2)]);
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    export_fig('./FIGS/ALLPCE_BBZ_S1s_0.png', '-dpng', '-r300')
    
    yyaxis right
    for k=1:size(Z_S1,2)
        semilogx(Rcofs(:,1), Z_S1(:,k), '-', 'Color', colos(k,:), 'LineWidth', 2)
    end
    ylabel('First Order Sobol Indices (sq. root)')
    ylim([-0.1 1.1])
%     xlim(minmax(Rcofs(:,1)'))
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_BBZ_S1s.png', '-dpng', '-r300')
    %%
    figure(1000)
    clf()
    set(gcf, 'Color', 'white')
    colos = DISTINGUISHABLE_COLORS(size(W_S1,2));
    ff = gobjects(length(ps),1);
    for j=1:length(ps)
        ff(j) = fill([Rcofs(:,1); Rcofs(end:-1:1,1)]*0, ...
            [WPI(:, j); WPI(end:-1:1, j+length(ps))]/2/pi*0, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
        legend(ff(j), sprintf('$%d^{th}-%d^{th}$ Percentiles', fix(ps(j)*100), fix((1-ps(j))*100)));
    end
    set(gca, 'XScale', 'log')
    aa = gobjects(size(W_S1,2), 1);
    for j=1:length(exp)
        aa(1) = semilogx(exp(j).AMP_avg*0, exp(j).FRE_avg*0, 'k-', 'LineWidth', 2);
    end
    legend(aa(1), 'Experimental Data')
    xlabel('Response Amplitude (m)')
    ylabel('Natural Frequency (Hz)')
    
    yyaxis right
    for k=1:size(W_S1,2)
        aa(1+k) = semilogx(Rcofs(:,1)*0, W_S1(:,k)*0, '-', 'Color', colos(k,:), 'LineWidth', 2);
        legend(aa(1+k), sprintf('$S_{%s}$', Pars{is{i}(k)}))
    end
    ll=legend([ff;aa]);
    set(ll, 'NumColumns', 3);
    set(ll, 'Position', [0.2 0.5 0.6325 0.1924])
    axis off
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_BBL_S1s.png', '-dpng', '-r300')
    
    %% S2 Plots
    % Frequency - S2
    figure((i-1)*10+2)
    clf()
    set(gcf, 'Color', 'white')
    colos = DISTINGUISHABLE_COLORS(size(W_S2,2));
    for j=1:length(ps)
        ff(j) = fill([Rcofs(:,1); Rcofs(end:-1:1,1)], ...
            [WPI(:, j); WPI(end:-1:1, j+length(ps))]/2/pi, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
        legend(ff(j), sprintf('$%d^{th}-%d^{th}$ Percentiles', fix(ps(j)*100), fix((1-ps(j))*100)));
    end
    set(gca, 'XScale', 'log')
    aa = gobjects(size(W_S2,2)+1, 1);
    for j=1:length(exp)
        aa(1)=semilogx(exp(j).AMP_avg, exp(j).FRE_avg, 'k-', 'LineWidth', 2);
    end
    ll=legend(aa(1), 'Experimental Data'); set(ll, 'Visible', 'off')
    xlabel('Response Amplitude (m)')
    ylabel('Natural Frequency (Hz)')
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    export_fig('./FIGS/ALLPCE_BBW_S2s_0.png', '-dpng', '-r300')
    
    yyaxis right
    for k=1:size(W_S2,2)
        aa(1+k) = semilogx(Rcofs(:,1), W_S2(:,k), '-', 'Color', colos(k,:), 'LineWidth', 2);
        legend(aa(1+k), sprintf('$S_{%s,%s}$', Pars{j1j2s(k,1)}, Pars{j1j2s(k,2)}))
    end
    ll=legend([ff;aa], 'Location', 'west');
    set(ll, 'Visible', 'off')
    ylabel('Second Order Sobol Indices (sq. root)')
    ylim([-0.1 1.1])
%     xlim(minmax(Rcofs(:,1)'))
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_BBW_S2s.png', '-dpng', '-r300')
    
    % Damping - S2
    figure((i-1)*10+6)
    clf()
    set(gcf, 'Color', 'white')
    for j=1:length(ps)
        fill([Rcofs(:,1); Rcofs(end:-1:1,1)], ...
            [ZPI(:, j); ZPI(end:-1:1, j+length(ps))]*100+zeta0*100, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
    end
    set(gca, 'XScale', 'log')
    for j=1:length(exp)
        semilogx(exp(j).AMP_avg, exp(j).DAM_avg*100, 'k-', 'LineWidth', 2);
    end
    xlabel('Response Amplitude (m)')
    ylabel('Damping Factor (\%)')
    yl = ylim;
    ylim([0 yl(2)]);
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    export_fig('./FIGS/ALLPCE_BBZ_S2s_0.png', '-dpng', '-r300')
    
    yyaxis right
    for k=1:size(W_S2,2)
        semilogx(Rcofs(:,1), Z_S2(:,k), '-', 'Color', colos(k,:), 'LineWidth', 2)
    end
    ylabel('Second Order Sobol Indices (sq. root)')
    ylim([-0.1 1.1])
%     xlim(minmax(Rcofs(:,1)'))
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_BBZ_S2s.png', '-dpng', '-r300')
    
    figure(1000)
    clf()
    set(gcf, 'Color', 'white')
    colos = DISTINGUISHABLE_COLORS(size(W_S2,2));
    for j=1:length(ps)
        ff(j) = fill([Rcofs(:,1); Rcofs(end:-1:1,1)]*0, ...
            [WPI(:, j); WPI(end:-1:1, j+length(ps))]/2/pi*0, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
        legend(ff(j), sprintf('$%d^{th}-%d^{th}$ Percentiles', fix(ps(j)*100), fix((1-ps(j))*100)));
    end
    set(gca, 'XScale', 'log')
    aa = gobjects(size(W_S2,2)+1, 1);
    for j=1:length(exp)
        aa(1)=semilogx(exp(j).AMP_avg*0, exp(j).FRE_avg*0, 'k-', 'LineWidth', 2);
    end
    legend(aa(1), 'Experimental Data')
    xlabel('Response Amplitude (m)')
    ylabel('Natural Frequency (Hz)')
    
    yyaxis right
    for k=1:size(W_S2,2)
        aa(1+k) = semilogx(Rcofs(:,1)*0, W_S2(:,k)*0, '-', 'Color', colos(k,:), 'LineWidth', 2);
        legend(aa(1+k), sprintf('$S_{%s,%s}$', Pars{j1j2s(k,1)}, Pars{j1j2s(k,2)}))
    end
    ll=legend([ff;aa], 'Location', 'west');
    ylabel('Second Order Sobol Indices (sq. root)')
    ylim([-0.1 1.1])
%     xlim(minmax(Rcofs(:,1)'))
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    set(gcf, 'Color', 'white')
    axis off
    set(ll, 'NumColumns', 3);
    set(ll, 'Position', [0.19 0.5 0.6325 0.1924])
    export_fig('./FIGS/ALLPCE_BBL_S2s.png', '-dpng', '-r300')
    
    %% ST1 Plots    
    % Frequency - ST1
    figure((i-1)*10+3)
    clf()
    set(gcf, 'Color', 'white')
    colos = DISTINGUISHABLE_COLORS(size(W_ST1,2));
    ff = gobjects(length(ps),1);
    for j=1:length(ps)
        ff(j) = fill([Rcofs(:,1); Rcofs(end:-1:1,1)], ...
            [WPI(:, j); WPI(end:-1:1, j+length(ps))]/2/pi, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
        legend(ff(j), sprintf('$%d^{th}-%d^{th}$ Percentiles', fix(ps(j)*100), fix((1-ps(j))*100)));
    end
    set(gca, 'XScale', 'log')
    aa = gobjects(size(W_ST1,2), 1);
    for j=1:length(exp)
        aa(1) = semilogx(exp(j).AMP_avg, exp(j).FRE_avg, 'k-', 'LineWidth', 2);
    end
    ll=legend(aa(1), 'Experimental Data'); set(ll, 'Visible', 'off')
    xlabel('Response Amplitude (m)')
    ylabel('Natural Frequency (Hz)')
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
%     export_fig('./FIGS/ALLPCE_BBW_STs_0.png', '-dpng', '-r300')
    print('./FIGS/ALLPCE_BBW_STs_0.png', '-dpng', '-r300')
    
    yyaxis right
    for k=1:size(W_ST1,2)
        aa(1+k) = semilogx(Rcofs(:,1), W_ST1(:,k), '-', 'Color', colos(k,:), 'LineWidth', 2);
        legend(aa(1+k), sprintf('$S_{T%s}$', Pars{is{i}(k)}))
    end
    ll=legend([ff;aa], 'Location', 'southwest');
    set(ll, 'Visible', 'off')
    ylabel('Total Sobol Indices (sq. root)')
    ylim([-0.1 1.1])
%     xlim(minmax(Rcofs(:,1)'))
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    set(gcf, 'Color', 'white')
%     export_fig('./FIGS/ALLPCE_BBW_STs.png', '-dpng', '-r300')
    print('./FIGS/ALLPCE_BBW_STs.png', '-dpng', '-r300')
    
    % Damping - ST1
    figure((i-1)*10+7)
    clf()
    set(gcf, 'Color', 'white')
    for j=1:length(ps)
        fill([Rcofs(:,1); Rcofs(end:-1:1,1)], ...
            [ZPI(:, j); ZPI(end:-1:1, j+length(ps))]*100+zeta0*100, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
    end
    set(gca, 'XScale', 'log')
    for j=1:length(exp)
        semilogx(exp(j).AMP_avg, exp(j).DAM_avg*100, 'k-', 'LineWidth', 2);
    end
    xlabel('Response Amplitude (m)')
    ylabel('Damping Factor (\%)')
    yl = ylim;
    ylim([0 yl(2)]);
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
%     export_fig('./FIGS/ALLPCE_BBZ_STs_0.png', '-dpng', '-r300')
    print('./FIGS/ALLPCE_BBZ_STs_0.png', '-dpng', '-r300')
    
    yyaxis right
    for k=1:size(Z_ST1,2)
        semilogx(Rcofs(:,1), Z_ST1(:,k), '-', 'Color', colos(k,:), 'LineWidth', 2)
    end
    ylabel('Total Sobol Indices (sq. root)')
    ylim([-0.1 1.1])
%     xlim(minmax(Rcofs(:,1)'))
    xlim([min(Rcofs(:,1)) min(10^(-4.5), max(Rcofs(:,1)))])
    set(gcf, 'Color', 'white')
%     export_fig('./FIGS/ALLPCE_BBZ_STs.png', '-dpng', '-r300')
    print('./FIGS/ALLPCE_BBZ_STs.png', '-dpng', '-r300')
    
    figure(1000)
    clf()
    set(gcf, 'Color', 'white')
    colos = DISTINGUISHABLE_COLORS(size(W_S1,2));
    ff = gobjects(length(ps),1);
    for j=1:length(ps)
        ff(j) = fill([Rcofs(:,1); Rcofs(end:-1:1,1)]*0, ...
            [WPI(:, j); WPI(end:-1:1, j+length(ps))]/2/pi*0, fcs(j,:), ...
            'EdgeColor', ecs(j,:), 'FaceAlpha', falph, 'LineStyle', lst{j}); hold on
        legend(ff(j), sprintf('$%d^{th}-%d^{th}$ Percentiles', fix(ps(j)*100), fix((1-ps(j))*100)));
    end
    set(gca, 'XScale', 'log')
    aa = gobjects(size(W_S1,2), 1);
    for j=1:length(exp)
        aa(1) = semilogx(exp(j).AMP_avg*0, exp(j).FRE_avg*0, 'k-', 'LineWidth', 2);
    end
    legend(aa(1), 'Experimental Data')
    xlabel('Response Amplitude (m)')
    ylabel('Natural Frequency (Hz)')
    
    yyaxis right
    for k=1:size(W_S1,2)
        aa(1+k) = semilogx(Rcofs(:,1)*0, W_S1(:,k)*0, '-', 'Color', colos(k,:), 'LineWidth', 2);
        legend(aa(1+k), sprintf('$S_{T%s}$', Pars{is{i}(k)}))
    end
    ll=legend([ff;aa]);
    set(ll, 'NumColumns', 3);
    set(ll, 'Position', [0.2 0.5 0.6325 0.1924])
    axis off
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_BBL_STs.png', '-dpng', '-r300')
    
    %% Modal Frequencies
    figure((i-1)*10+4)
    clf()
    set(gcf, 'Color', 'white')
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
%     set(gca, 'yscale', 'log')
    ll = legend(bb, 'Location', 'best');
    set(ll, 'Visible', 'off')
    xlabel('Bending Mode ID')
    ylabel('Sobol Indices (sq. root)')
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_MFREQ_S1s.eps', '-depsc')
    
    figure(1000)
    clf()
    set(gcf, 'Color', 'white')
    colos = DISTINGUISHABLE_COLORS(size([Wstat_S1 Wstat_S2],2));
%     bb = bar([Wstat_S1(mdis, :) Wstat_S2(mdis, :)]*nan);
    bb = gobjects(size(Wstat_S1,2)+size(Wstat_S2,2),1);
    for k=1:size(Wstat_S1,2)
        bb(k) = fill(zeros(10,1), zeros(10,1), colos(k,:)); hold on
%         bb(k).FaceColor = colos(k,:);
        legend(bb(k), sprintf('$S_{%s}$', Pars{is{i}(k)}));
    end
    for k=size(Wstat_S1,2)+(1:size(Wstat_S2,2))
        bb(k) = fill(zeros(10,1), zeros(10,1), colos(k,:));
%         bb(k).FaceColor = colos(k,:);
        legend(bb(k), sprintf('$S_{%s, %s}$', Pars{j1j2s(k-size(Wstat_S1,2),1)}, Pars{j1j2s(k-size(Wstat_S1,2),2)}));
    end
%     set(gca, 'yscale', 'log')
    axis off
    ll = legend(bb);
    set(ll, 'NumColumns', 5)
    set(ll, 'Position', [0.1 0.5 0.8 0.24])
    xlabel('Bending Mode ID')
    ylabel('Sobol Indices (sq. root)')
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_MFREQL_S1s.eps', '-depsc')
    %%
    figure((i-1)*10+8)
    clf()
    set(gcf, 'Color', 'white')
    for im=1:length(mdis)
        subplot(length(mdis), 1, im)
        histogram(rmoutliers(Psixs*Wstatcofs(mdis(im),:)'/2/pi), 'Normalization', 'pdf'); hold on
        [f, x] = ksdensity(rmoutliers(Psixs*Wstatcofs(mdis(im),:)'/2/pi));
        plot(x, f, '-', 'LineWidth', 2); hold on
        ylabel('PDF')
        title(sprintf('Mode %d', im))
%         ylim([0 0.75])
        plot(Wsex(im)/2/pi*[1 1], ylim, 'k-', 'LineWidth', 2)
    end
    xlabel('Natural Frequency (Hz)')
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_MFREQsims.eps', '-depsc')
    
    %%
    figure((i-1)*10+9)
    clf()
    set(gcf, 'Color', 'white')
    colos = DISTINGUISHABLE_COLORS(size(Wstat_ST1,2));
    bb = bar(Wstat_ST1(mdis, :));
    for k=1:size(Wstat_ST1,2)
        bb(k).FaceColor = colos(k,:);
        legend(bb(k), sprintf('$S_{T%s}$', Pars{is{i}(k)}));
    end
%     set(gca, 'yscale', 'log')
    legend(bb, 'Location', 'best')
    xlabel('Bending Mode ID')
    ylabel('Total Sobol Indices (sq. root)')
    set(gcf, 'Color', 'white')
    export_fig('./FIGS/ALLPCE_MFREQ_STs.eps', '-depsc')
end
