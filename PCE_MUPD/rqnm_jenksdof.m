clc
clear all
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/QUADRATURE/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

%% Model Parameters
m = 1.0;
k = 4.0;
c = 2*0.02*sqrt(k/m);


fnl = @(x, kt,muN) (kt*x).*(abs(kt*x)<muN)+sign(kt*x).*muN.*(abs(kt*x)>=muN);
Dnl = @(x, kt,muN) (4*muN.*(abs(x)-muN/kt)).*(abs(kt*x)>=muN);

%% RQNM along the backbone
lam = @(x, kt,muN) x.*(k*x+fnl(x, kt,muN))./(m*x.^2);
ztx = @(x, kt,muN) (Dnl(x, kt,muN)+pi*c*x.^2.*lam(x, kt,muN))./(2*pi*x.^2.*lam(x, kt,muN).^1.5);

%% Gaussian Parameters
kt = makedist('normal', 5, 0.1);
muN = makedist('normal', 0.1, 0.01);

% %% Gamma Distribution Parameters -> Doesn't Work
% kt = makedist('gam', 2500, 1/500);
% muN = makedist('gam', 100, 1/1000);
% 
% %% Gamma Germs
% quadptfun = @(n) LAGWT(n);
% polfun = @(n, x) PLAGU(n, x);
% tformfun = @(Xdatr, kt, muN) Xdatr.*[kt.b muN.b];
% tpref = 'LAGU';

% %% Uniform germs
% quadptfun = @(n) LGWT(n, -1, 1);
% polfun = @(n, x) PLEGE(n, x);
% tformfun = @(Xdatr, kt, muN) icdf('normal', (Xdatr+[1 1])/2, [kt.mu muN.mu], [kt.sigma muN.sigma]);
% % tformfun = @(Xdatr, kt, muN) [icdf(kt, (Xdatr(:,1)+1)/2) icdf(muN, (Xdatr(:,2)+1)/2)];
% % tformfun = @(Xdatr, kt, muN) [icdf(kt, cdf('Uniform', Xdatr(:,1), -1, 1)) icdf(muN, cdf('Uniform', Xdatr(:,2), -1, 1))];
% tpref = 'LAGP';

%% Gaussian Germs
quadptfun = @(n) GPHWT(n);
polfun = @(n, x) PHERM(n, x);
tformfun = @(Xdatr, kt, muN) Xdatr.*[kt.sigma muN.sigma] + [kt.mu muN.mu];
tpref = 'HERM';

%% Create samples from quadrature points
Nq = 10;
[xi, wi] = quadptfun(Nq);

[Xi, Yi] = meshgrid(xi, xi);
Wi = reshape(wi.*wi', [], 1);
Xdatr = [Xi(:) Yi(:)];

%% Transform to normal distribution
Xdat = tformfun(Xdatr, kt, muN);

N = Nq^2;

%% Simulations
Nx = 200;
xs = logspace(-4, 2, Nx);

Wsdat = zeros(N, Nx);
Zsdat = zeros(N, Nx);
for i=1:N
    Wsdat(i, :) = sqrt(lam(xs, Xdat(i,1), Xdat(i,2)));
    Zsdat(i, :) = ztx(xs, Xdat(i,1), Xdat(i,2));
end

% "mean model"
Wsmm = sqrt(lam(xs, mean(kt), mean(muN)));
Zsmm = ztx(xs, mean(kt), mean(muN));
%% PCE Basis evaluation
Np = 5;

IJs = [repmat((0:Np)', Np+1,1) kron((0:Np)', ones(Np+1,1))];
IJs = IJs(sum(IJs,2)<=Np,:);

P = size(IJs,1);
Psi = zeros(Nq^2, P);
Integs = zeros(P, 1);  % Not really needed, just used for checking quadrature
for n=1:P
    [p1, i1] = polfun(IJs(n,1), Xdatr(:, 1));
    [p2, i2] = polfun(IJs(n,2), Xdatr(:, 2));
    
    Psi(:, n) = p1.*p2;
    Integs(n) = i1*i2;
end
% Check if (Wi'*(Psi.^2)) is equal to (Integs')
% It will be scaled by 2*2=4 for LEGPOL, i.e., the area of the reference
% square ((-1,-1) to (1,1)).
Integs = (Wi'*(Psi.^2))';

ix = find(IJs(2:end,2)==0)+1;  % Only x
iy = find(IJs(2:end,1)==0)+1;  % Only x
ixy = setdiff((2:P), [ix; iy]);  % x & y

%% Estimate PCE coefficients at each amplitude level through quadrature integration
W_As = zeros(P, Nx);
Z_As = zeros(P, Nx);

% Interpolation Errors
W_errs = zeros(Nx, 1);
Z_errs = zeros(Nx, 1);

% Total variance of the PC models
W_DPC = zeros(Nx,1);
Z_DPC = zeros(Nx,1);

% First order sobol indices
W_S1 = zeros(Nx, 2);  
Z_S1 = zeros(Nx, 2);

% Second order sobol indices
W_S2 = zeros(Nx, 1);
Z_S2 = zeros(Nx, 1);
for ia=1:Nx
    W_As(:, ia) = (Wi'*(Wsdat(:, ia).*Psi))./(Integs');
    Z_As(:, ia) = (Wi'*(Zsdat(:, ia).*Psi))./(Integs');
    
    % Interpolation Errors
    W_errs(ia) = rms(Wsdat(:, ia)-Psi*W_As(:, ia));
    Z_errs(ia) = rms(Zsdat(:, ia)-Psi*Z_As(:, ia));
    
    W_DPC(ia) = (W_As(2:end, ia).^2)'*Integs(2:end);
    Z_DPC(ia) = (Z_As(2:end, ia).^2)'*Integs(2:end);
    
    W_S1(ia, :) = [((W_As(ix, ia).^2)'*Integs(ix))/W_DPC(ia) ((W_As(iy, ia).^2)'*Integs(iy))/W_DPC(ia)];
    Z_S1(ia, :) = [((Z_As(ix, ia).^2)'*Integs(ix))/Z_DPC(ia) ((Z_As(iy, ia).^2)'*Integs(iy))/Z_DPC(ia)];
    
    W_S2(ia) = ((W_As(ixy, ia).^2)'*Integs(ixy))/W_DPC(ia);
    Z_S2(ia) = ((Z_As(ixy, ia).^2)'*Integs(ixy))/Z_DPC(ia);
end

%% Plots
figure(10)
clf()
set(gcf, 'color','white')

subplot(2,1, 1)
bb1 = semilogx(xs, Wsdat, '.', 'Color', [1 1 1]*0.8); hold on
% semilogx(xs, mean(Wsdat), 'k-', 'LineWidth', 1); hold on;
bb2 = semilogx(xs, Wsmm, 'c-', 'LineWidth', 2); hold on
bb3 = semilogx(xs, W_As(1, :), 'k-', 'LineWidth', 1); hold on;
% semilogx(xs, mean(Wsdat)+sqrt(var(Wsdat)).*[-1; 1]*3, '-', 'Color', [1 1 1]*0.6);
bb4 = semilogx(xs, W_As(1,:)+sqrt(W_DPC').*[-1; 1]*3, '-.', 'Color', [1 1 1]*0);
ylabel('Natural Frequency')

yyaxis right; 
semilogx(xs, abs(W_S1(:,1)), 'b-', xs, abs(W_S1(:,2)), 'r-', xs, abs(W_S2), 'm-'); 
set(gca, 'yscale', 'linear')
ylabel('$SU_i$')
ylim([-0.1 1.1])

legend([bb1(1) bb2(1) bb3(1) bb4(1)], 'Data', 'Mean Model', 'PC Mean', 'PC 6-$\sigma$ Variance Bound', 'Location', 'west')
if strcmp(tpref, 'HERM')    
    title('Hermite Polynomial')
elseif strcmp(tpref, 'LAGP')
    title('Lagrange Polynomial')
else
    title('Laguerre Polynomial')
end


subplot(2,1, 2)
semilogx(xs, Zsdat, '.', 'Color', [1 1 1]*0.8); hold on
% semilogx(xs, mean(Zsdat), 'k-'); hold on;
semilogx(xs, Zsmm, 'c-', 'LineWidth', 2); hold on
semilogx(xs, Z_As(1,:), 'k-', 'LineWidth', 1); hold on;
% semilogx( xs, mean(Zsdat)+sqrt(var(Zsdat)).*[-1; 1]*3, '-', 'Color', [1 1 1]*0.6);
semilogx(xs, Z_As(1,:)+sqrt(Z_DPC').*[-1; 1]*3, '-.', 'Color', [1 1 1]*0);

ylabel('Damping Factor')
yyaxis right; 
aa=semilogx(xs, abs(Z_S1(:,1)), 'b-', xs, abs(Z_S1(:,2)), 'r-', xs, abs(Z_S2), 'm-');
set(gca, 'yscale', 'linear')
ylabel('$SU_i$')
ylim([-0.1 1.1])
xlabel('Resonant Amplitude')

if strcmp(tpref, 'HERM')
    legend(aa, '$SU_{k_t}$', '$SU_{\mu N}$', '$SU_{k_t,\mu N}$', 'Location', 'west')
end

% export_fig(sprintf('./FIGS/%s_respind.png', tpref), '-dpng')

figure(20)
clf()
set(gcf, 'color','white')
subplot(2,1, 1);
loglog(xs, W_As(1,:), 'k-'); hold on
loglog(xs, W_errs, 'k-.')
ylabel('Natural Frequency')
% ylim([1e-15 1e2])

if strcmp(tpref, 'HERM') 
    title('Hermite Polynomial')
else
    title('Lagrange Polynomial')
end

subplot(2,1, 2);
loglog(xs, Z_As(1,:), 'k-'); hold on
loglog(xs, Z_errs, 'k-.')
ylabel('Damping Factor')
xlabel('Resonant Amplitude')
% ylim([1e-15 1e2])

% export_fig(sprintf('./FIGS/%s_resperr.png', tpref), '-dpng')