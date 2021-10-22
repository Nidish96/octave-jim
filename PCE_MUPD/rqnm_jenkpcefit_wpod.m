clc
clear all
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

fnl = @(x, kt, muN) (kt*x).*(abs(kt*x)<muN)+sign(kt*x).*muN.*(abs(kt*x)>=muN);
jnl = @(x, kt, muN) (kt).*(abs(kt*x)<muN);
Dnl = @(x, kt, muN) (4*muN.*(abs(x)-muN/kt)).*(abs(kt*x)>=muN);

%% RQNM along the backbone
lam = @(x, klin, clin, kt,muN) x.*(klin*x+fnl(x, kt,muN))./(m*x.^2);
ztx = @(x, klin, clin, kt,muN) (Dnl(x, kt,muN).*(lam(x, klin,clin,kt,muN))+pi*clin*x.^2.*lam(x, klin,clin,kt,muN))./(2*m*pi*x.^2.*lam(x, klin,clin,kt,muN).^1.5);

%% Gaussian Parameters
klin = makedist('normal', Wmin_exp^2, Wmin_exp^2/10);
clin = makedist('normal', 2*Zt_exp*Wmin_exp, 2*Zt_exp*Wmin_exp/10);

kt = makedist('normal', (Wstat_exp^2-klin.mu), (Wstat_exp^2-klin.mu)/10);
muN = makedist('normal', kt.mu*Qs(end/2), kt.mu*Qs(end));

%% Mean Model
Wmm = sqrt(lam(Qs, klin.mu, clin.mu, kt.mu, muN.mu));
Zmm = ztx(Qs, klin.mu, clin.mu, kt.mu, muN.mu);

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
tformfun = @(Xdatr, klin, clin, kt, muN) abs(Xdatr.*[klin.sigma clin.sigma kt.sigma muN.sigma] + [klin.mu clin.mu kt.mu muN.mu]);
tpref = 'HERM';

% %% Uniform germs
% quadptfun = @(n) LGWT(n, -1, 1);
% polfun = @(n, x) PLEGE(n, x);
% tformfun = @(Xdatr, klin, clin, kt, muN) icdf('normal', (Xdatr+[1 1 1 1])/2, [klin.mu clin.mu kt.mu muN.mu], [klin.sigma clin.sigma kt.sigma muN.sigma]);
% % tformfun = @(Xdatr, kt, muN) [icdf(kt, (Xdatr(:,1)+1)/2) icdf(muN, (Xdatr(:,2)+1)/2)];
% % tformfun = @(Xdatr, kt, muN) [icdf(kt, cdf('Uniform', Xdatr(:,1), -1, 1)) icdf(muN, cdf('Uniform', Xdatr(:,2), -1, 1))];
% tpref = 'LAGP';

%% Create samples from quadrature points
Nq = 10;
[xi, wi] = quadptfun(Nq);

[Xi1, Xi2, Xi3, Xi4] = ndgrid(xi);
[Wi1, Wi2, Wi3, Wi4] = ndgrid(wi);

Xis = [Xi1(:) Xi2(:) Xi3(:) Xi4(:)];
Xds = tformfun(Xis, klin, clin, kt, muN);
Wis = Wi1(:).*Wi2(:).*Wi3(:).*Wi4(:);

clear Xi1 Xi2 Xi3 Xi4 Wi1 Wi2 Wi3 Wi4

% Get rid of design where the quadrature weight is less than threshold
thresh = 1e-15;
epi = find(Wis/max(Wis)<thresh);
Xis(epi, :) = [];
Xds(epi, :) = [];
Wis(epi, :) = [];
Nqp = size(Xis,1);

%% Setup PCE
Np = 5;  % Number of basis functions
[I1, I2, I3, I4] = ndgrid(0:Np);
Is = [I1(:) I2(:) I3(:) I4(:)];
Is = Is(sum(Is, 2)<=Np, :);

clear I1 I2 I3 I4

% Construct Polynomials at Quadrature Points
Npol = size(Is, 1);  % Number of polynomial bases;

Psiqp = ones(Nqp, Npol);
Integs = ones(Npol, 1);
for i=1:4
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
    Wldats(is) = sqrt(Xds(is,1)+jnl(0, Xds(is, 3), Xds(is, 4)));
    Zldats(is) = Xds(is, 3)/(2*Wldats(is));
    Wdats(is, :) = sqrt(lam(Qs, Xds(is, 1), Xds(is, 2), Xds(is, 3), Xds(is, 4)));
    Zdats(is, :) = ztx(Qs, Xds(is, 1), Xds(is, 2), Xds(is, 3), Xds(is, 4));
    
%     fprintf('%d/%d\n', is, Nqp)
end

% %% Calculate POD of Data
% thresh = 1e-15;
% [U, S, V] = svd([Wdats Zdats], 'econ');
% Nr = find(diag(S)/S(1,1)>thresh, 1, 'last' );  % Number of reduced outputs
% Lr = V(:, 1:Nr);
% % Datsr = [Wdats Zdats]*Lr;
% Datsr = U(:, 1:Nr)*S(1:Nr, 1:Nr);

%% Use Raw Data
Nr = Nx*2;
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

%% Initialize RLS
P_0 = eye(size(D_As, 1))*1e-3;
Dcs_0 = D_As;

rls_lamb = 1;  % forgetfulness

Npop = 100;
xs = [random(klin, Npop, 1) random(clin, Npop, 1) random(kt, Npop, 1) random(muN, Npop, 1)];
xs(1,:) = 0;

NITs = 5;

Objs = cell(NITs, 1);
Sols = cell(NITs, 1); 

for it=1:NITs
    %%  Conduct Optimization with gamultiobj
    costfun = @(x) PCECOSTFUN(x, eye(4), Dcs_0, Lr(1:Nx,:), Lr(Nx+1:end,:), [Wse(:), Zse(:)], Is, repmat({polfun}, 1, 4));

    % costfun(kron((0:2)'/2, ones(1, 4)));
    figure(1);

    opt = optimoptions('gamultiobj', 'PopulationSize', Npop, ...
        'UseVectorized', true, 'UseParallel', true, ...
        'InitialPopulationMatrix', xs, ...
        'PlotFcn', {'gaplotpareto', 'gaplotparetodistance', 'gaplotrankhist', 'gaplotspread'});

    [x, Objs{it}, ~, op] = gamultiobj(@(x) costfun(x), ...
        4, [], [], [], [], zeros(1, 4), [], [], opt);
    xd = tformfun(x, klin, clin, kt, muN);

    %% Check Against original model simulations
    Wtrue = zeros(size(x,1), Nx);
    Ztrue = zeros(size(x,1), Nx);

    [Dpce, Psi] = EVALPCE(x, Dcs_0', Is, repmat({polfun}, 1, 4));
    Dpce = Dpce';
    for is=1:size(x,1)
        Wtrue(is, :) = sqrt(lam(Qs, xd(is, 1), xd(is, 2), xd(is, 3), xd(is, 4)));
        Ztrue(is, :) = ztx(Qs, xd(is, 1), xd(is, 2), xd(is, 3), xd(is, 4));
    end
    Dtrue = [Wtrue Ztrue]*Lr;

    %% Add the new points to the previous using recursive least-squares (uniformly weighted)
    errmax = 0;
    for is=1:length(x)
        P_t = (P_0 - (P_0*(Psi(is,:)'*Psi(is,:))*P_0)/(rls_lamb+Psi(is,:)*P_0*Psi(is,:)'))/rls_lamb;

        Dcs_t = Dcs_0 + P_t*Psi(is, :)'*(Dtrue(is,:)-Psi(is,:)*Dcs_0);

        Dis_0 = Psi(is, :)*Dcs_0;
        Dis_t = Psi(is, :)*Dcs_t;

        errmax = max(errmax, rms(Dis_0-Dis_t));

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
    % Updated Model
    Dpce1 = EVALPCE(x, Dcs_0', Is, repmat({polfun}, 1, 4));
    
    %% Save
    Sols{it} = [Wtrue Ztrue Dpce*Lr Dpce1*Lr];

    %% Plot
%     figure(10)
%     clf()
%     subplot(2,1, 1)
%     % semilogx(Qs, Dpce*Lr(1:Nx, :)'/2/pi, '.-', Qs, Wse/2/pi, 'k-')
%     semilogx(Qs, Dtrue*Lr(1:Nx, :)'/2/pi, '.-', Qs, Wse/2/pi, 'k-')
%     subplot(2,1, 2)
%     % semilogx(Qs, Dpce*Lr(Nx+(1:Nx), :)', '.-', Qs, Zse, 'k-')
%     semilogx(Qs, Dtrue*Lr(Nx+(1:Nx), :)', '.-', Qs, Zse, 'k-')

    %% Check Convergence
    fprintf('%d. Error = %e\n', it, errmax);
end

%% Plots



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
%         zres(i) = rms((Zs(i, :)'-WZexp(:,2)));  %./WZexp(:,2));
        zres(i) = rms((log10(Zs(i, :)')-log10(WZexp(:,2)))./log10(WZexp(:,2)));
    end
    
    wzres = [wres, zres];
end
