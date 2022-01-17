clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/QUADRATURE/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

mdi = 1;
%% Load Experimental Data
exp(1) = load(sprintf('./DATA/Mode%d_Low.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');
exp(2) = load(sprintf('./DATA/Mode%d_Med.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');
exp(3) = load(sprintf('./DATA/Mode%d_High.mat', mdi), 'AMP_avg', ...
              'FRE_avg', 'DAM_avg');

fav = robustfit([exp(1).AMP_avg; exp(2).AMP_avg; exp(3).AMP_avg], ...
                [exp(1).FRE_avg; exp(2).FRE_avg; exp(3).FRE_avg]);
zav = robustfit([exp(1).AMP_avg; exp(2).AMP_avg; exp(3).AMP_avg], ...
                [exp(1).DAM_avg; exp(2).DAM_avg; exp(3).DAM_avg]);
Wstat_exp = 2*pi*fav(1);
Zt_exp    = zav(1);

Wmin_exp = 2*pi*mean(cellfun(@(c) c(end), {exp.FRE_avg}));

[Qe, si] = sort(cell2mat({exp.AMP_avg}'));
We = 2*pi*cell2mat({exp.FRE_avg}');  We = We(si);
Ze = cell2mat({exp.DAM_avg}');  Ze = Ze(si);
Qe = Qe./We.^2;
si = find(isfinite(Qe));
Qe = Qe(si);
We = We(si);
Ze = Ze(si);

%% Smooth Data
Nx = 100;
Qs = logspace(log10(min(Qe))+0.1, log10(max(Qe))-0.1, Nx);

Wsm = smoothdata(We, 'sgolay');
Zsm = smoothdata(Ze, 'sgolay');

Wse = interp1(Qe, Wsm, Qs);
Zse = interp1(Qe, Zsm, Qs);

figure(1)
clf()
subplot(2,1,1)
semilogx(Qe, We/2/pi, '.'); hold on
semilogx(Qe, Wsm/2/pi, 'k-', 'LineWidth', 2)
% plot(Qs, Wse/2/pi, 'ro')
ylabel('Natural Frequency (Hz)')
subplot(2,1,2)
semilogx(Qe, Ze, '.'); hold on
semilogx(Qe, Zsm, 'k-', 'LineWidth', 2)
% plot(Qs, Zse, 'ro')
ylabel('Damping Factor')
xlabel('Response Amplitude (m)')

%% Model Parameters
m = 1.0;

phimax = @(fs, kt, chi, bt) abs(fs).*(1+abs(bt))./(abs(kt).*(abs(bt)+(chi+1)./(chi+2 )));
rfun = @(fs, kt, chi, bt) abs(fs).*(chi+1)./(phimax(abs(fs),abs(kt),chi,abs(bt)).^(chi+2).*(abs(bt)+(chi+1)./(chi+2)));
sfun = @(fs, kt, chi, bt) bt.*rfun(abs(fs),abs(kt),chi,abs(bt)).*phimax(abs(fs),abs(kt),chi,abs(bt)).^(chi+1)./(chi+1);

fnl = @(x, fs, kt, chi, bt) (abs(kt).*x - rfun(fs, kt, chi, bt).*x.^(chi+2)./((chi+1).*(chi+2))).*(abs(x)<phimax(fs,kt,chi,bt))+abs(fs).*sign(x).*(abs(x)>=phimax(fs,kt,chi,bt));
jnl = @(x, fs, kt, chi, bt) (abs(kt) - rfun(fs, kt, chi, bt).*x.^(chi+1)./(chi+1)).*(abs(x)<phimax(fs,kt,chi,bt));

% phimax = @(fs, kt, chi, bt) fs.*(1+bt)./(kt.*(bt+(chi+1)./(chi+2 )));
% rfun = @(fs, kt, chi, bt) fs.*(chi+1)./(phimax(fs,kt,chi,bt).^(chi+2).*(bt+(chi+1)./(chi+2)));
% sfun = @(fs, kt, chi, bt) bt.*rfun(fs,kt,chi,bt).*phimax(fs,kt,chi,bt).^(chi+1)./(chi+1);
% 
% fnl = @(x, fs, kt, chi, bt) (kt.*x - rfun(fs, kt, chi, bt).*x.^(chi+2)./((chi+1).*(chi+2))).*(abs(x)<phimax(fs,kt,chi,bt))+fs.*sign(x).*(abs(x)>=phimax(fs,kt,chi,bt));
% jnl = @(x, fs, kt, chi, bt) (kt - rfun(fs, kt, chi, bt).*x.^(chi+1)./(chi+1)).*(abs(x)<phimax(fs,kt,chi,bt));

Dnl = @(x, fs, kt, chi, bt) (4*rfun(fs,kt,chi,bt).*x.^(chi+3)./((chi+2).*(chi+3))).*(abs(x)<phimax(fs,kt,chi,bt)) + (4*rfun(fs,kt,chi,bt).*phimax(fs,kt,chi,bt).^(chi+2).*(x./(chi+2)-phimax(fs,kt,chi,bt)./(chi+3))+4*sfun(fs,kt,chi,bt).*phimax(fs,kt,chi,bt).*(x-phimax(fs,kt,chi,bt))).*(abs(x)>=phimax(fs,kt,chi,bt));

%% RQNM along the backbone
lam = @(x, pars) x.*(abs(pars(1))*x+fnl(x, pars(3),pars(4),pars(5),pars(6)))./(m*x.^2);
ztx = @(x, pars) (Dnl(x, pars(3),pars(4),pars(5),pars(6)).*sqrt(lam(x, pars))+pi*abs(pars(2))*x.^2.*lam(x, pars))./(2*m*pi*x.^2.*lam(x, pars).^1.5);


% lam = @(x, klin, clin, fs,kt,chi,bt) x.*(klin*x+fnl(x, fs,kt,chi,bt))./(m*x.^2);
% ztx = @(x, klin, clin, fs,kt,chi,bt) (Dnl(x, fs,kt,chi,bt).*sqrt(lam(x, klin,clin,fs,kt,chi,bt))+pi*clin*x.^2.*lam(x, klin,clin,fs,kt,chi,bt))./(2*m*pi*x.^2.*lam(x, klin,clin,fs,kt,chi,bt).^1.5);

%% Gaussian Parameters
% klin = makedist('normal', Wstat_exp^2, Wstat_exp^2/10);
% clin = makedist('normal', 2*Zt_exp*Wstat_exp, 2*Zt_exp*Wstat_exp/10);

% kt = makedist('normal', 0, (Wstat_exp^2-Wmin_exp^2));
% muN = makedist('normal', 0, (Wstat_exp^2-Wmin_exp^2)/100);

klin = makedist('normal', Wmin_exp^2, Wmin_exp^2/20);
clin = makedist('normal', 2*Zt_exp*Wmin_exp, 2*Zt_exp*Wmin_exp/20);

kt = makedist('normal', (Wstat_exp^2-klin.mu), (Wstat_exp^2-klin.mu)/20);
fs = makedist('normal', kt.mu*Qs(end/2), kt.mu*Qs(end));
chi = makedist('normal', 0, 1/3);
bt = makedist('normal', 0, 1e-4);

%% Mean Model
Wmm = sqrt(lam(Qs, [klin.mu, clin.mu, fs.mu, kt.mu, chi.mu, bt.mu]));
Zmm = ztx(Qs, [klin.mu, clin.mu, fs.mu, kt.mu, chi.mu, bt.mu]);

figure(1)
subplot(2,1,1); 
hold on
semilogx(Qs, Wmm/2/pi, '.-')
subplot(2,1,2); 
hold on
semilogx(Qs, Zmm, '.-')

%% Gaussian Germs
quadptfun = @(n) GPHWT(n);
polfun = @(n, x) PHERM(n, x);
tformfun = @(Xdatr, klin, clin, fs, kt, chi, bt) (Xdatr.*[klin.sigma clin.sigma fs.sigma kt.sigma chi.sigma bt.sigma] + [klin.mu clin.mu fs.mu kt.mu chi.mu bt.mu]);
tpref = 'HERM';

%% Uniform germs
% quadptfun = @(n) LGWT(n, -1, 1);
% polfun = @(n, x) PLEGE(n, x);
% tformfun = @(Xdatr, klin, clin, fs, kt, chi, bt) icdf('normal', (Xdatr+ones(1,6))/2, [klin.mu clin.mu fs.mu kt.mu chi.mu bt.mu], [klin.sigma clin.sigma fs.sigma kt.sigma chi.sigma bt.sigma]);
% % tformfun = @(Xdatr, kt, muN) [icdf(kt, (Xdatr(:,1)+1)/2) icdf(muN, (Xdatr(:,2)+1)/2)];
% % tformfun = @(Xdatr, kt, muN) [icdf(kt, cdf('Uniform', Xdatr(:,1), -1, 1)) icdf(muN, cdf('Uniform', Xdatr(:,2), -1, 1))];
% tpref = 'LAGP';

%% Create samples from quadrature points
Nq = 5;
[xi, wi] = quadptfun(Nq);

[Xi1, Xi2, Xi3, Xi4, Xi5, Xi6] = ndgrid(xi);
[Wi1, Wi2, Wi3, Wi4, Wi5, Wi6] = ndgrid(wi);

Xis = [Xi1(:) Xi2(:) Xi3(:) Xi4(:) Xi5(:) Xi6(:)];
Xds = tformfun(Xis, klin, clin, fs, kt, chi, bt);
Wis = Wi1(:).*Wi2(:).*Wi3(:).*Wi4(:).*Wi5(:).*Wi6(:);

clear Xi1 Xi2 Xi3 Xi4 Xi5 Xi6 Wi1 Wi2 Wi3 Wi4 Wi5 Wi6

% Get rid of design where the quadrature weight is less than threshold
% thresh = 1e-15;
% epi = find(Wis/max(Wis)<thresh);
% Xis(epi, :) = [];
% Xds(epi, :) = [];
% Wis(epi, :) = [];
Nqp = size(Xis,1);

disp(Nqp)

%% Setup PCE
Np = 5;  % Number of basis functions
[I1, I2, I3, I4, I5, I6] = ndgrid(0:Np);
Is = [I1(:) I2(:) I3(:) I4(:) I5(:) I6(:)];
Is = Is(sum(Is, 2)<=Np, :);

clear I1 I2 I3 I4 I5 I6

% Construct Polynomials at Quadrature Points
Npol = size(Is, 1);  % Number of polynomial bases;

Psiqp = ones(Nqp, Npol);
Integs = ones(Npol, 1);
for i=1:6
    [pfn, pint] = polfun(Is(:, i), Xis(:, i));
    
    Psiqp = Psiqp.*pfn;
    Integs = Integs.*pint;
end
% Integs = (Wi'*(Psi.^2))';  % Numerically determine the integrals

%% Conduct Simulations
Wldats = zeros(Nqp, 1);
Zldats = zeros(Nqp, 1);
Wdats = zeros(Nqp, Nx);
Zdats = zeros(Nqp, Nx);
for is=1:Nqp
    Wldats(is) = sqrt(Xds(is,1)+jnl(0, Xds(is, 3), Xds(is, 4), Xds(is, 5), Xds(is, 6)));
    Zldats(is) = Xds(is, 3)/(2*Wldats(is));
    Wdats(is, :) = sqrt(lam(Qs, Xds(is, :)));
    Zdats(is, :) = ztx(Qs, Xds(is, :));
    
    fprintf('%d/%d\n', is, Nqp)
end

% %% Calculate POD of Data
% thresh = 1e-15;
% [U, S, V] = svd([Wdats Zdats], 'econ');
% Nr = find(diag(S)/S(1,1)>thresh, 1, 'last' );  % Number of reduced outputs
% Lr = V(:, 1:Nr);
% % Datsr = [Wdats Zdats]*Lr;
% Datsr = U(:, 1:Nr)*S(1:Nr, 1:Nr);

%%
Nr = 2*Nx;
Datsr = [Wdats Zdats];
Lr = eye(Nr);

%% Estimate PCE Coefficients
D_As = zeros(Npol, Nr);

% Quadrature
for ia=1:Nr
    D_As(:, ia) = (Wis'*(Datsr(:, ia).*Psiqp))./(Integs');
end

% % Regression
% D_As = (Psiqp'*Psiqp)\Psiqp'*Datsr;

%%
Dpce = EVALPCE(zeros(1, 6), D_As', Is, repmat({polfun}, 1, 6));
Dpce = Dpce';

figure(1)
subplot(2,1,1); 
hold on
semilogx(Qs, Dpce*Lr(1:Nx, :)'/2/pi, 'o-')
subplot(2,1,2); 
hold on
semilogx(Qs, Dpce*Lr(Nx+(1:Nx),:)', 'o-')

%% Optimize True Model
pars0 = [klin.mu, clin.mu, fs.mu, kt.mu, chi.mu, bt.mu];

r2lfun = @(pars) [log10(pars(1:4)) pars(5:6)];
l2rfun = @(lpars) [10.^lpars(1:4) lpars(5:6)];
lbs = [-1, -1, -1, -1, -1, 0];
ubs = [inf, inf, inf, inf, 1, inf];

% r2lfun = @(pars) pars;
% l2rfun = @(lpars) lpars;
% lbs = [0, 0, 0, 0, -1, 0];
% ubs = [inf, inf, inf, inf, 1, inf];

resfun = @(pars) [rms(sqrt(lam(Qs, l2rfun(pars)))./Wse-1), ...
    rms(log10(ztx(Qs, l2rfun(pars)))./log10(Zse)-1)];

% %%
Npop = 100;
opt = optimoptions('gamultiobj', 'PopulationSize', Npop, ...
    'UseVectorized', false, 'UseParallel', true, ...
    'PlotFcn', {'gaplotpareto', 'gaplotparetodistance', 'gaplotrankhist', 'gaplotspread'});
    
rng(1);
inpop = abs(randn(Npop, 6).*[klin.sigma clin.sigma fs.sigma kt.sigma chi.sigma bt.sigma]+[klin.mu, clin.mu, fs.mu, kt.mu, chi.mu, bt.mu]);
for i=1:Npop
    inpop(i, :) = r2lfun(inpop(i,:));
end
opt.InitialPopulationMatrix = inpop;

[parsol, objsol] = gamultiobj(@(x) resfun(x), ...
    6, [], [], [], [], lbs, ubs, [], opt);

trSolz = zeros(Nx, size(parsol,1)*2);
for i=1:size(parsol,1)
    trSolz(:, i+[0 size(parsol,1)]) = [sqrt(lam(Qs, l2rfun(parsol(i,:))))' ztx(Qs, l2rfun(parsol(i,:)))'];
end

%% Initialize RLS
P_0 = eye(size(D_As, 1))*1e-1;  %1e-3
Dcs_0 = D_As;

rls_lamb = 0.9;  % forgetfulness

rng(1)

Npop = 100;
% xds = [random(klin, Npop, 1) random(clin, Npop, 1) random(fs, Npop, 1) random(kt, Npop, 1) random(chi, Npop, 1) random(bt, Npop, 1)];
% xds(1:2,:) = 0;
% xds(2, 3) = -0.01;

% xs = (xds-[klin.mu clin.mu fs.mu kt.mu chi.mu bt.mu])./[klin.sigma clin.sigma fs.sigma kt.sigma chi.sigma bt.sigma];
% xs(1, :) = 0;

xs = zeros(1, 6);

NITs = 20;  % Five Iterations
errs = zeros(NITs,1);
Objs = cell(NITs, 1);
Sols = cell(NITs, 1); 

opt = optimoptions('gamultiobj', 'PopulationSize', Npop, ...
    'UseVectorized', true, 'UseParallel', true, ...
    'InitialPopulationMatrix', xs, ...
    'PlotFcn', {'gaplotpareto', 'gaplotparetodistance', 'gaplotrankhist', 'gaplotspread'});
    
for it=1:NITs
    %%  Conduct Optimization with gamultiobj
    costfun = @(x) PCECOSTFUN(x, eye(6), Dcs_0, Lr(1:Nx,:), Lr(Nx+1:end,:), [Wse(:), Zse(:)], Is, repmat({polfun}, 1, 6));

    % costfun(kron((0:2)'/2, ones(1, 4)));
    figure(1);

    opt.InitialPopulationMatrix = xs;
    [xs, Objs{it}, ~, op] = gamultiobj(@(x) costfun(x), ...
        6, [], [], [], [], zeros(1, 6), [], [], opt);
    xd = tformfun(xs, klin, clin, fs, kt, chi, bt);

    %% Check Against original model simulations
    Wtrue = zeros(size(xs,1), Nx);
    Ztrue = zeros(size(xs,1), Nx);

    [Dpce, Psi] = EVALPCE(xs, Dcs_0', Is, repmat({polfun}, 1, 6));
    Dpce = Dpce';
    for is=1:size(xs,1)
        Wtrue(is, :) = sqrt(lam(Qs, xd(is, :)));
        Ztrue(is, :) = ztx(Qs, xd(is, :));
    end
    Dtrue = [Wtrue Ztrue]*Lr;

    %% Add the new points to the previous using recursive least-squares (uniformly weighted)
    errmax = 0;
    for is=1:length(xs)
        P_t = (P_0 - (P_0*(Psi(is,:)'*Psi(is,:))*P_0)/(rls_lamb+Psi(is,:)*P_0*Psi(is,:)'))/rls_lamb;

        Dcs_t = Dcs_0 + P_t*Psi(is, :)'*(Dtrue(is,:)-Psi(is,:)*Dcs_0);

        Dis_0 = Psi(is, :)*Dcs_0;
        Dis_t = Psi(is, :)*Dcs_t;

        errmax = max(errmax, rms((Dis_0-Dis_t)./Dis_t));

    %     figure(1); clf(); 
    %     subplot(2,1,1);
    %     semilogx(Qs, Wse, 'k-'); hold on
    %     semilogx(Qs, Wtrue(is,:), 'k-'); hold on
    %     semilogx(Qs, Dis_0*Lr(1:Nx, :)', 'bo'); hold on
    %     semilogx(Qs, Dis_t*Lr(1:Nx, :)', 'r-', 'LineWidth', 2); hold on
    %     title(sprintf('Point %d', is))
    %     subplot(2,1,2); 
    %     semilogx(Qs, Zse, 'k-'); hold on
    %     semilogx(Qs, Ztrue(is,:), 'k-'); hold on
    %     semilogx(Qs, Dis_0*Lr(Nx+1:end, :)', 'bo'); hold on
    %     semilogx(Qs, Dis_t*Lr(Nx+1:end, :)', 'r-', 'LineWidth', 2); hold on    
    %     pause(0.5)

        P_0 = P_t;
        Dcs_0 = Dcs_t;
    end
    errs(it) = errmax;
    
    % Updated Model
    Dpce1 = EVALPCE(xs, Dcs_0', Is, repmat({polfun}, 1, 6));
    Dpce1 = Dpce1';
    
    %% Save
    Sols{it} = [Wtrue Ztrue Dpce*Lr Dpce1*Lr];

    %% Plot
    figure(10)
    clf()
    subplot(2,1, 1)
    % semilogx(Qs, Dpce*Lr(1:Nx, :)'/2/pi, '.-', Qs, Wse/2/pi, 'k-')
    semilogx(Qs, Dtrue*Lr(1:Nx, :)'/2/pi, '.-', Qs, Wse/2/pi, 'k-')
    subplot(2,1, 2)
    % semilogx(Qs, Dpce*Lr(Nx+(1:Nx), :)', '.-', Qs, Zse, 'k-')
    semilogx(Qs, Dtrue*Lr(Nx+(1:Nx), :)', '.-', Qs, Zse, 'k-')

    %% Error
    fprintf('%d. Error = %e\n', it, errmax);
    
end

%% Plots
[~, im] = min(errs);  % Minimum error iteration

figure(1)
clf()
subplot(2,1, 1)
semilogx(Qs, Sols{im}(:, 1:Nx)/2/pi, '-'); hold on
semilogx(Qs, Wse/2/pi, 'k-', 'LineWidth', 3);
xlim([2e-7 1e-5])
ylabel('Natural Frequency (Hz)')
subplot(2,1, 2)
semilogx(Qs, Sols{im}(:, Nx+(1:Nx))*100, '-'); hold on
semilogx(Qs, Zse*100, 'k-', 'LineWidth', 3)
xlim([2e-7 1e-5])
xlabel('Amplitude (m)')
ylabel('Damping Factor (%)')

set(gcf, 'Color', 'white')
% export_fig('./FIGS/converged_bbs.eps', '-depsc')

%%
cols = DISTINGUISHABLE_COLORS(im);

figure(110)
clf()
set(gcf, 'Color', 'white')
aa = gobjects(im+1, 1);
[~, si] = sort(objsol(:,1));
aa(1) = plot(objsol(si, 1), objsol(si, 2), 'k-', 'LineWidth', 2); hold on
legend(aa(1), 'True NDS')
for it=1:im
    [~, si] = sort(Objs{it}(:,1));
    
    aa(1+it) = plot(Objs{it}(si,1), Objs{it}(si, 2), 'h', ...
        'Color', cols(it,:), 'MarkerFaceColor', cols(it,:)); hold on
    legend(aa(it+1), sprintf('It. %d', it))
    legend(aa(1:(it+1)), 'Location', 'eastoutside');
    
    xlabel('Rel-RMS Frequency Deviation')
    ylabel('Rel-RMS Log-Damping Deviation')
    
%     export_fig(sprintf('./FIGS/PARANIM/paranim_%02d.png', it), '-dpng')
end
% export_fig('./FIGS/Paretofronts.eps', '-depsc')

%%
figure(12)
clf()
semilogy(errs(1:im), 'ko-', 'MarkerFaceColor', 'w')
xlabel('Iteration Index');
ylabel('Max. Relative Surrogate Model Deviation')

set(gcf, 'Color', 'white')
% export_fig('./FIGS/convergence.eps', '-depsc')

%% PCE Cost Function
function [wzres] = PCECOSTFUN(x, Lx, Dcofs, Lrw, Lrz, WZexp, Is, polfuns)    
    if size(x, 2)~=size(Is,2)
        x = x';
    end
    x = x*Lx';
    
    % PCE Bases
    Psi = ones(size(x, 1), size(Is, 1));
    for j=1:size(Is, 2)
        psi = polfuns{j}(Is(:,j), x(:, j));
        Psi = Psi.*psi;
    end
    Ds = Psi*Dcofs;
    Ws = Ds*Lrw';
    Zs = Ds*Lrz';
    
    wres = zeros(size(x,1),1);
    zres = zeros(size(x,1),1);
    for i=1:size(x,1)
        wres(i) = rms((Ws(i, :)'-WZexp(:,1))./WZexp(:,1));
%         zres(i) = rms((Zs(i, :)'-WZexp(:,2))); % ./WZexp(:,2));
        zres(i) = rms((log10(Zs(i, :)')-log10(WZexp(:,2)))./log10(WZexp(:,2)));
    end
    
    wzres = [wres, zres];
end