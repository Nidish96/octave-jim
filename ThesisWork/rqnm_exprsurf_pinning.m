clc
% clear all

addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/GRAPH/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/CONTACTMODELS/')
addpath('../ROUTINES/QUASISTATIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

model = 'BRB_Thesis';
E = 1.9231e11;
nu = 0.3;

% model = 'BRB';
% E = 2e11;
% nu = 0.3;

top   = 'R05B_After';
bot   = 'R05A_After';

load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'K', 'M', 'R', 'L', 'Fv', 'TFM');

%% Load Mesh
Nds = dlmread(sprintf('../MODELS/%s/Nodes.dat', model));
Quad = dlmread(sprintf('../MODELS/%s/Elements.dat', model));

Nq = 2;
MESH = MESH2D(Nds, 3, [], Quad, Nq);

%% Prepare Contact Model Parameters
MESH = MESH.SETQUAD(1);
Aels = full(sum(MESH.Tm));  % Element Areas
Aint = sum(Aels);
Aels = kron(Aels(:), ones(Nq^2,1));
MESH = MESH.SETQUAD(Nq);

load(sprintf('./MATFILES/%s_R1_AspPDEs.mat', top), 'PS_sds');
R1top = load(sprintf('./MATFILES/%s_R1_AspPDEs.mat', top), 'BilinPlaneQPs', 'LLX0s_sd', 'CRAD', 'NASPS', 'PS_sds');
R2top = load(sprintf('./MATFILES/%s_R2_AspPDEs.mat', top), 'BilinPlaneQPs', 'LLX0s_sd', 'CRAD', 'NASPS');

R1bot = load(sprintf('./MATFILES/%s_R1_AspPDEs.mat', bot), 'BilinPlaneQPs', 'LLX0s_sd', 'CRAD', 'NASPS');
R2bot = load(sprintf('./MATFILES/%s_R2_AspPDEs.mat', bot), 'BilinPlaneQPs', 'LLX0s_sd', 'CRAD', 'NASPS');

%% Gap Function ( with Nqp = 2^2 per element )
gap1 = R1top.BilinPlaneQPs(:,1)-R1bot.BilinPlaneQPs(:,1);
gap2 = R2top.BilinPlaneQPs(:,1)-R2bot.BilinPlaneQPs(:,1);
gap1 = gap1-max(gap1);
gap2 = gap2-max(gap2);

gap = (gap1+gap2)/2;
gapsd = sqrt((R1top.BilinPlaneQPs(:,2).^2+R1bot.BilinPlaneQPs(:,2).^2+R2top.BilinPlaneQPs(:,2).^2+R2bot.BilinPlaneQPs(:,2).^2)/4);  % Gap Standard Deviation
%% Nasps, Lambda and z0
Nasps1 = R1top.NASPS(:,1)+R1bot.NASPS(:,2);
Nasps2 = R2top.NASPS(:,1)+R2bot.NASPS(:,2);
Nasps = (Nasps1+Nasps2)/2;

lam1 = (R1top.NASPS(:,1)+R1bot.NASPS(:,1))./(R1top.NASPS(:,1)./R1top.LLX0s_sd(:,1)+R1bot.NASPS(:,1)./R1bot.LLX0s_sd(:,1));
lam2 = (R2top.NASPS(:,1)+R2bot.NASPS(:,1))./(R2top.NASPS(:,1)./R2top.LLX0s_sd(:,1)+R2bot.NASPS(:,1)./R2bot.LLX0s_sd(:,1));
lam = (lam1+lam2)/2;

z01 = -log(0.01)./lam1;
z02 = -log(0.01)./lam2;
z0 = (z01+z02)/2;

% Collect
Nasps = kron(Nasps, ones(Nq^2,1));
lam   = kron(lam, ones(Nq^2,1));
% z0    = kron(z0, ones(Nq^2,1));
% z0 = -log(0.1)./lam;
z0 = log(Nasps)./lam;

%% Curvature Radii
R1 = (R1top.CRAD(:,1).*R1top.NASPS(:,1)+R1bot.CRAD(:,1).*R1bot.NASPS(:,1))./(R1top.NASPS(:,1)+R1bot.NASPS(:,2));
R2 = (R2top.CRAD(:,1).*R2top.NASPS(:,1)+R2bot.CRAD(:,1).*R2bot.NASPS(:,1))./(R2top.NASPS(:,1)+R2bot.NASPS(:,2));

Rad = (R1+R2)/2;
Rad = kron(Rad, ones(Nq^2,1));

%% Prestress
Prestress = mean([12002 12075 12670]);

%% Coefficient of Friction (AISI 304N SS)
% s = 186e6; H = 655e6;
% s = 190e6;  H = 545e6;
% s = 2.00;  H = 1.0;
s = 0.85;  H = 1.0;
mu = s/H*ones(MESH.Ne*MESH.Nq^2,1);

% mu=0.85;
% Z. Guo and M. R. W. Brake, “Uncertainty Quantification of Friction Model Parameters for Joint Mechanics Applications,” 2018 STLE Annual Meeting, Minneapolis, MN, May, 2018. 

%% Contact Model
Estar = E/(1-nu^2)/2;
Gstar = E/((1+nu)*(2-nu))/2;

% Estar = Estar/2;
% Gstar = Gstar/2;

cn = Estar*sqrt(pi*Rad./lam.^3).*Nasps./Aels;
ct = 4*Gstar*sqrt(pi*Rad./lam).*Nasps./Aels;

lamo = lam;  muo = mu; gapo = gap;
cno = cn.*exp(-lam.*z0); cto = ct.*exp(-lam.*z0);

%% Detect Bolt Nodes
% Adjacency Matrix
[~, Adj] = NODEELADJ(MESH.Nds, MESH.Quad(:, 2:end));

% Find bolt nodes as those that're not the convexhull
nis = find(diag(Adj)==2);
nis = setdiff(nis, nis(convhull(MESH.Nds(nis,:))));
% nis = nis(abs(MESH.Nds(nis, 1))<0.05 & abs(MESH.Nds(nis, 2))<0.005);

% Cluster the nodes into 3 bolt holes
rng(1)
bis = kmeans(MESH.Nds(nis,:), 3);
nis = [nis(bis==1)';
       nis(bis==2)';
       nis(bis==3)'];

% Transformation Matrix
Tb = zeros(3, MESH.Nn);
for i=1:3
    Tb(i, nis(i,:)) = 1/length(nis);
end

% figure(2); 
% clf()
% plot(graph(Adj-diag(diag(Adj))), 'XData', MESH.Nds(:,1), 'YData', MESH.Nds(:,2)); hold on
% for i=1:3
%     plot(MESH.Nds(nis(i,:), 1), MESH.Nds(nis(i,:), 2), '*')
% end

%% Create Object
GM = MDOFGEN(M, K, zeros(size(M)), L);

GM = GM.SETNLFUN(2+5, ...
    kron(MESH.Qm, eye(3))*L(1:(MESH.Nn*3), :), ...
    @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lamo, muo, gapo, varargin{:}), ...
    L(1:(MESH.Nn*3), :)'*kron(MESH.Tm, eye(3)));
GM = GM.SETNLFUN(1+3, ...
    kron(Tb, eye(2))*L(setdiff(1:(MESH.Nn*3), 3:3:(MESH.Nn*3)), :), ...
    @(t, u, ud) BPINNING(t, u, ud, pi/4*Estar*25.4e-3, 0.433e-3/2));

%% Linearized Stiffness
knlin = lamo.*Prestress/Aint;
ktlin = cto./cno*Prestress/Aint;

K0 = zeros(size(L,1));
K0(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = MESH.Tm*diag(ktlin)*MESH.Qm;
K0(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = MESH.Tm*diag(ktlin)*MESH.Qm;
K0(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = MESH.Tm*diag(knlin)*MESH.Qm;
K0 = L'*K0*L;

%% Prestress Analysis
fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');

% load('rqnmsolmm1.mat', 'Ustat');
% U0 = Ustat;
% [Ustat, ~, eflag, ~, Jstat] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
% if eflag<=0  % No convergence
    U0 = (K+K0)\(Fv*Prestress);
    U0 = L*U0;  U0(3:3:MESH.Nn*3) = MESH.Qm\gapo;
    U0 = L\U0;

    GM.NLTs(1).func = @(t, u, varargin) EXPROUGHFRICT(t, u, cto*0, cno, lamo, muo, gapo, varargin{:});
    [U0, ~, eflag] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
    if eflag<=0
        load('rqnmsolmm2.mat', 'Ustat');
        U0 = Ustat;
        U0 = L*U0;  U0(3:3:MESH.Nn*3) = MESH.Qm\gapo;
        U0 = L\U0;

        [U0, ~, eflag] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
        if eflag<=0
            error('No Prestress Convergence')
        end
    end
    GM.NLTs(1).func = @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lamo, muo, gapo, varargin{:});
    [Ustat, ~, eflag, ~, Jstat] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
% end

% %% Plot Static Tractions
% Tstat = GM.NLTs(1).func(0, GM.NLTs(1).L*Ustat);
% 
% figure(3)
% clf()
% for i=1:3
%     subplot(3,1,i)
%     MESH.SHOWFIELD2D(Tstat(i:3:end))
%     xx=colorbar('southoutside');
%     xlabel(xx, sprintf('Traction %d', i))
%     axis equal
% end
% colormap(jet)

%% Linearized Analysis
[Vstat, Wstat] = eigs(Jstat, M, 10, 'SM');
[Wstat, si] = sort(sqrt(diag(Wstat)));
Vstat = Vstat(:, si);
Vstat = Vstat./sqrt(diag(Vstat'*M*Vstat))';

%% March
mds = [1 3 5];

% mdi = 1;  % Mode of Interest
% AMIN = -0.5;  AMAX = 2.5;

% mdi = 2;  % Mode of Interest
% AMIN = 0;  AMAX = 3;

% mdi = 3;  % Mode of Interest
% AMIN = -2;  AMAX = 3;

Na = 10;
As = logspace(AMIN, AMAX, Na)*9.81/Wstat(mds(mdi))^2;
As = [-As(end:-1:1) As]';
Eflags = zeros(1, 2*Na);
UlC = zeros(GM.Ndofs+1, 2*Na);
dUdalC = zeros(GM.Ndofs+1, 2*Na);

tic
ul0 = [Ustat+Vstat(:,mds(mdi))*As(Na+1); Wstat(mds(mdi))^2];
opts = struct('reletol', 1e-6, 'Display', true, 'lsrch', 0, 'ITMAX', 50);
fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
for ia=1:Na    
    [UlC(:, Na+ia), ~, Eflags(Na+ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(As(Na+ia))], 1, Fv*Prestress, Ustat), ul0, fopts);
    eflag = Eflags(Na+ia);
    if eflag<=0
        error('No Convergence')
    end
    
    dRdA = [zeros(GM.Ndofs,1);-2*UlC(end, Na+ia)*As(Na+ia)]; 
    dUdalC(:, Na+ia) = -dRdUl\dRdA;
    if ia<Na
        ul0 = UlC(:, Na+ia);
        ul0(1:end-1) = Ustat+(ul0(1:end-1)-Ustat)*As(Na+ia+1)/As(Na+ia);
    end
    
    fprintf('%d, ', ia);
end
fprintf('\n');

ul0 = [Ustat+Vstat(:,mds(mdi))*As(Na+1-1); Wstat(mds(mdi))^2];
for ia=1:Na
    [UlC(:, Na+1-ia), ~, Eflags(Na+1-ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(-As(Na+1-ia))], 1, Fv*Prestress, Ustat), ul0, fopts);
    eflag = Eflags(Na+1-ia);
    if eflag<=0
        error('No Convergence')
    end
    
    dRdA = [zeros(GM.Ndofs,1);-2*UlC(end, Na+1-ia)*As(Na+1-ia)];
    dUdalC(:, Na+1-ia) = -dRdUl\dRdA;
    if ia<Na
        ul0 = UlC(:, Na+1-ia);
        ul0(1:end-1) = Ustat+(ul0(1:end-1)-Ustat)*As(Na+1-ia-1)/As(Na+1-ia);
    end
    
    fprintf('%d, ', ia);
end
toc

%% Post-Processing with zero harmonics
Ln = reshape([UlC(end,:); dUdalC(end,:)], 2*size(UlC,2), 1);
Un = reshape([permute(UlC(1:end-1,:), [3, 2, 1]); 
    permute(dUdalC(1:end-1,:), [3, 2, 1])], [], GM.Ndofs);  % (2Npt,Ndofs)

Nq = 100;  %% TESTING
% Qs = (10^Amax)*linspace(0.01, 1, Nq)';
Qs = logspace(QMIN(mdi), QMAX(mdi), Nq)';

Nt = 2^7;
t = linspace(0, 2*pi, Nt+1)'; t(end) = [];

qw20 = [0; Wstat(1)^2];
qw2s = zeros(3, Nq);

opt = optimset('Jacobian', 'on', 'Display', 'off');

qi = 1;
dUt = zeros(Nt, Nq, GM.Ndofs);
Ut = zeros(Nt, Nq, GM.Ndofs);
Udot = zeros(Nt, Nq, GM.Ndofs);
for qi=[Nq:-1:1 2:Nq]
    qw2s(1:end-1, qi) = fsolve(@(qw2) RQHBRESFUN([qw2; Qs(qi)], @(t, q) HERMINTERP(As, Ln, q(:)).*q(:), ...
        @(t, q) HERMINTERP(As, Ln, q(:))+interp1(As, Ln(2:2:end), q).*q(:), Nt), qw20, opt);
    qw2s(end, qi) = Qs(qi);
    qw20 = qw2s(1:end-1, qi);
    
    % Mode Shapes
    qt = qw2s(1,qi)+qw2s(3,qi)*cos(t);
    [~, Nint, dNint] = HERMINTERP(As, Ln, qt);
    
    Ut(:, qi, :) = Nint*Un;
    dUt(:, qi, :) = dNint*Un;
    Udot(:, qi, :) = (dNint*Un).*(-sqrt(qw2s(2,qi))*qw2s(3,qi)*sin(t));
    
    fprintf('%d\n', qi)
end
Lams = qw2s(end-1, :);
qdot = -sin(t).*(sqrt(Lams(:)).*Qs(:))';
qddot = -cos(t).*(Lams(:).*Qs(:))';

% Uh = reshape(GETFOURIERCOEFF(1, reshape(Ut, Nt, Nq*GM.Ndofs)), 2, Nq, GM.Ndofs);
% Phi = (squeeze(Uh(1,:,:)-1j*Uh(2,:,:))./Qs)';

dUh = squeeze(reshape(GETFOURIERCOEFF(0, reshape(dUt, Nt, Nq*GM.Ndofs)), 1, Nq, GM.Ndofs))';
Phi = dUh;

Phi = Phi./sqrt(diag(Phi'*GM.M*Phi)');

%% Damping
tol = 1e-6;
Nits = 2;  % Maximum marching iterations
tic
Zts = zeros(Nq, 1);
tic
parfor qi=1:Nq
%     Zts(qi) = GETFOURIERCOEFF(0, sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M +...
%         squeeze(Udot(:, qi, :))*GM.C +...
%         squeeze(Ut(:, qi, :))*GM.K + ...
%         GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)), tol, Nits)),2))/(Qs(qi)^2*Lams(qi)^1.5);
    Zts(qi) = GETFOURIERCOEFF(0, sum(squeeze(Udot(:, qi, :)).*(squeeze(Udot(:, qi, :))*GM.C +...
        GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)), tol, Nits)),2))/(Qs(qi)^2*Lams(qi)^1.5);
    fprintf('%d\n', qi)
end
toc

save(sprintf('./MEANMODELRESP/MMP_MODE%d.mat', mdi), 'Qs', 'Lams', 'Zts', 'Phi', 'As', 'UlC', 'dUdalC', 'Ustat', 'Wstat')
%% Plotting
Rxs = abs(R(3,:)*Phi);
Rx = R(3,:)*(UlC(1:end-1,:)-Ustat)./As';

% load('./MATFILES/SaveFile_2021-Feb-19_High.mat', 'AMP_avg', 'FRE_avg', 'DAM_avg');
load(sprintf('./MATFILES/Mode%d_High.mat', mdi), 'AMP_avg', 'FRE_avg', 'DAM_avg');

figure((mdi-1)*100+1)
clf()
set(gcf, 'Color', 'white')
semilogx((Qs.*Rxs').*Lams(:)/9.81, sqrt(Lams)/2/pi, '-', 'LineWidth', 2); hold on
semilogx(AMP_avg, FRE_avg, 'k-', 'LineWidth', 2)

semilogx(abs(As.*Rx').*UlC(end,:)'/9.81, sqrt(UlC(end,:))/2/pi, 'ko', 'MarkerFaceColor', 'k')

legend('Model Prediction', 'Experimental Data', 'Location', 'Best')
xlabel('Response Amplitude (g)')
ylabel('Natural Frequency (Hz)')

xlim(minmax(((Qs.*Rxs').*Lams(:)/9.81)'))

% export_fig('./FIGS/BBWmm.png', '-dpng')

figure((mdi-1)*100+2)
clf()
set(gcf, 'Color', 'white')
semilogx((Qs.*Rxs').*Lams(:)/9.81, Zts*100, '-', 'LineWidth', 2); hold on
semilogx(AMP_avg, DAM_avg*100, 'k-', 'LineWidth', 2)

legend('Model Prediction', 'Experimental Data', 'Location', 'Best')
xlabel('Response Amplitude (g)')
ylabel('Damping Factor (%)')
xlim(minmax(((Qs.*Rxs').*Lams(:)/9.81)'))

% export_fig('./FIGS/BBZmm.png', '-dpng')