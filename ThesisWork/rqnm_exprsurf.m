clc
% clear all

addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
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

fnl = @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lamo, muo, gapo, varargin{:});

%% Create Object
GM = MDOFGEN(M, K, zeros(size(M)), L);

% Kt = [1e12; 1e12; 0];
% kn = 1e12;
% mu = 0.25;
% gap = 0;
% 
% fnl = @(t, u, varargin) ELDRYFRICT(t, u, Kt, kn, mu, gap, varargin{:});
GM = GM.SETNLFUN(2+5, ...
    kron(MESH.Qm, eye(3))*L(1:(MESH.Nn*3), :), ...
    fnl, ...
    L(1:(MESH.Nn*3), :)'*kron(MESH.Tm, eye(3)));

% Multiple Nonlinearities form - Do not use, VERY slow
% for iq=1:MESH.Nq^2*MESH.Ne
%     is = (iq-1)*3+(1:3);
%     
%     GM = GM.SETNLFUN(2+5, ...
%         kron(MESH.Qm(iq, :), eye(3))*L(1:(MESH.Nn*3), :), ...
%         @(t, u, varargin) EXPROUGHFRICT(t, u, cto(iq), cno(iq), lamo(iq), muo(iq), gapo(iq), varargin{:}), ...
%         L(1:(MESH.Nn*3), :)'*kron(MESH.Tm(:, iq), eye(3)));
%     
%     fprintf('%d\n', iq);
% end

%% Linearized Stiffness
knlin = lamo.*Prestress/Aint;
ktlin = cto./cno*Prestress/Aint;

K0 = zeros(size(L,1));
K0(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = MESH.Tm*diag(ktlin)*MESH.Qm;
K0(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = MESH.Tm*diag(ktlin)*MESH.Qm;
K0(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = MESH.Tm*diag(knlin)*MESH.Qm;
K0 = L'*K0*L;

%% Prestress Analysis
% load('statguess.mat', 'U0');
load('statguess1.mat', 'U0');
% load('newstat.mat', 'Ustat');  U0 = Ustat;

% U0 = (MESH.Qm*L(3:3:MESH.Nn*3,:))\gapo;
% U0 = (K+K0)\(Fv*Prestress);
% U0 = L*U0;  U0(3:3:MESH.Nn*3) = MESH.Qm\gapo;
% U0 = L\U0;

% fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% GM.NLTs.func = @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lamo, muo, gapo, varargin{:});
% [Ustat, ~, eflag, ~, Jstat] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);

%% Prestress Analysis
fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');

% load('rqnmsolmm1.mat', 'Ustat');
% U0 = Ustat;
% [Ustat, ~, eflag, ~, Jstat] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
% if eflag<=0  % No convergence
    U0 = (K+K0)\(Fv*Prestress);
    U0 = L*U0;  U0(3:3:MESH.Nn*3) = MESH.Qm\gapo;
    U0 = L\U0;

    GM.NLTs.func = @(t, u, varargin) EXPROUGHFRICT(t, u, cto*0, cno, lamo, muo, gapo, varargin{:});
    [U0, ~, eflag] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
    if eflag<=0
        load('rqnmsolmm1.mat', 'Ustat');
        U0 = Ustat;
        U0 = L*U0;  U0(3:3:MESH.Nn*3) = MESH.Qm\gapo;
        U0 = L\U0;

        [U0, ~, eflag] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
        if eflag<=0
            error('No Prestress Convergence')
        end
    end
    GM.NLTs.func = @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lamo, muo, gapo, varargin{:});
    [Ustat, ~, eflag, ~, Jstat] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
% end

%% Plot Static Tractions
Tstat = GM.NLTs.func(0, GM.NLTs.L*Ustat);

figure(2)
clf()
for i=1:3
    subplot(3,1,i)
    MESH.SHOWFIELD2D(Tstat(i:3:end))
    xx=colorbar('southoutside');
    xlabel(xx, sprintf('Traction %d', i))
    axis equal
end
colormap(jet)

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
        
%         [~, dRdUl, dRdlA] = GM.RQMRESFUN([UlC(:, Na+ia); log10(As(Na+ia))], 1, Fv*Prestress, Ustat);
%         duldlA = -dRdUl\dRdlA;
%         ul0 = ul0+duldlA*(log10(As(Na+ia+1))-log10(As(Na+ia)));
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
%         ul0 = [Ustat-(UlC(1:end-1, Na+1+ia)-Ustat); UlC(end-1,Na+1+ia)];
    end
    
    fprintf('%d, ', ia);
end
toc

% save('./rqnmsolmm1.mat', 'UlC', 'As', 'dUdalC', 'Ustat')

% %% GPU Computation
% 
% % fnl = @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lamo, muo, gapo, varargin{:});
% 
% uxs = gpuArray(GM.NLTs.L(1:3:end,:)*Ustat);
% uys = gpuArray(GM.NLTs.L(2:3:end,:)*Ustat);
% uns = gpuArray(GM.NLTs.L(3:3:end,:)*Ustat);
% ctg = gpuArray(cto);
% cng = gpuArray(cno);
% lamg = gpuArray(lamo);
% mug = gpuArray(muo);
% gapg = gpuArray(gapo);
% tg = gpuArray(zeros(size(uxs)));
% 
% % uxs = GM.NLTs.L(1:3:end,:)*Ustat;
% % uys = GM.NLTs.L(2:3:end,:)*Ustat;
% % uns = GM.NLTs.L(3:3:end,:)*Ustat;
% % ctg = cto;
% % cng = cno;
% % lamg = lamo;
% % mug = muo;
% % gapg = gapo;
% 
% tic
% Ts = arrayfun(@EXPROUGHFRICT_sca, tg, uxs, uys, uns, ctg, cng, lamg, mug, gapg);
% % us = arrayfun(@(u) u, uxs);
% toc
% % 
% % us = GM.NLTs.L*Ustat;
% % tic
% % Ts = sum(EXPROUGHFRICT(0, us, ctg, cng, lamg, mug, gapg));
% % toc

%%
OpDs = R(3,:)*(UlC(1:end-1,:)-Ustat);
OpAs = (R(3,:)*(UlC(1:end-1,:)-Ustat)).*(UlC(end,:))/9.81;
OpWs = sqrt(UlC(end,:));

figure(2)
clf()
semilogx(abs(OpAs(1:end)), OpWs(1:end)/2/pi, '.-')

figure(3)
clf()
plot(As, As.*UlC(end,:)', '.-')
bbu = As(Na+1:end)
bbf = As(Na+1:end).*UlC(end,Na+1:end)';
bbfm = As(Na:-1:1).*UlC(end,Na:-1:1)';
% Red & Blue are Masing Predictions
plot(bbu, bbf, 'ko-', ...
    bbu(end)-2*bbu, bbf(end)-2*bbf, 'b.-', ...
    -bbu(end)+2*bbu, -bbf(end)+2*bbf, 'r.-', ...
    -bbu, bbfm, '*-')

%% Post Processing (Hermite Interpolation + SHB)
Ln = reshape([UlC(end,:); dUdalC(end,:)], 2*size(UlC,2), 1);

QMIN = [-5.5 -6 -8.75];  QMAX = [-2.625 -3.16 -3.77];

Nq = 100;  %% TESTING
% Qs = (10^Amax)*linspace(0.01, 1, Nq)';
Qs = logspace(QMIN(mdi), QMAX(mdi), Nq)';

Nt = 2^7;
t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
qt = cos(t).*Qs';
[Lt, Nint, dNint] = HERMINTERP(As, Ln, qt(:));
Lt = reshape(Lt, Nt, Nq);

%% Natural Frequency
Lams = sqrt(sum((GETFOURIERCOEFF(1, Lt.*qt)./Qs').^2))';
qdot = -sin(t).*(sqrt(Lams).*Qs)';
qddot = -cos(t).*(Lams.*Qs)';

%% Mode Shape
Un = reshape([permute(UlC(1:end-1,:), [3, 2, 1]); 
    permute(dUdalC(1:end-1,:), [3, 2, 1])], [], GM.Ndofs);  % (2Npt,Ndofs)
Ut = reshape(Nint*Un, Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs)
Udot = reshape((dNint*Un).*qdot(:), Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs)
Uddot = reshape((dNint*Un).*qddot(:), Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs)

Uh = reshape(GETFOURIERCOEFF(1, reshape(Ut, Nt, Nq*GM.Ndofs)), 2, Nq, GM.Ndofs);
Phi = (squeeze(Uh(1,:,:)-1j*Uh(2,:,:))./Qs)';

%% Damping
tol = 1e-6;
Nits = 2;  % Maximum marching iterations
tic
Zts = zeros(Nq, 1);
tic
parfor (qi=1:Nq,8)
% for qi=1:Nq
%     size(sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M + squeeze(Udot(:, qi, :))*GM.C + squeeze(Ut(:, qi, :))*GM.K + GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)))),2))
    Zts(qi) = GETFOURIERCOEFF(0, sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M +...
        squeeze(Udot(:, qi, :))*GM.C +...
        squeeze(Ut(:, qi, :))*GM.K + ...
        GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)), tol, Nits)),2))/(Qs(qi)^2*Lams(qi)^1.5);
    fprintf('%d\n', qi)
end
toc

save(sprintf('./MEANMODELRESP/MM_MODE%d.mat', mdi), 'Qs', 'Lams', 'Zts', 'Phi', 'As', 'UlC', 'dUdalC', 'Ustat', 'Wstat')

%% Plotting
Rxs = abs(R(3,:)*Phi);
Rx = R(3,:)*(UlC(1:end-1,:)-Ustat)./As';

% load('./MATFILES/SaveFile_2021-Feb-19_High.mat', 'AMP_avg', 'FRE_avg', 'DAM_avg');
load(sprintf('./MATFILES/Mode%d_High.mat', mdi), 'AMP_avg', 'FRE_avg', 'DAM_avg');

figure(400)
clf()
set(gcf, 'Color', 'white')
semilogx((Qs.*Rxs').*Lams/9.81, sqrt(Lams)/2/pi, '-', 'LineWidth', 2); hold on
% semilogx(abs(As.*Rx').*UlC(end,:)'/9.81, sqrt(UlC(end,:))/2/pi, 'ko', 'MarkerFaceColor', 'k')
semilogx(AMP_avg, FRE_avg, 'k-', 'LineWidth', 2)

% xlim([10^-0.5 10^2.5])

legend('Model Prediction', 'Experimental Data', 'Location', 'Best')
xlabel('Response Amplitude (g)')
ylabel('Natural Frequency (Hz)')

% export_fig('./FIGS/BBWmm.png', '-dpng')

figure(500)
clf()
set(gcf, 'Color', 'white')
semilogx((Qs.*Rxs').*Lams/9.81, Zts*100, '-', 'LineWidth', 2); hold on
semilogx(AMP_avg, DAM_avg*100, 'k-', 'LineWidth', 2)

% xlim([10^-0.5 10^2.5])

legend('Model Prediction', 'Experimental Data', 'Location', 'Best')
xlabel('Response Amplitude (g)')
ylabel('Damping Factor (%)')
% ylim([-0.5 2])

% export_fig('./FIGS/BBZmm.png', '-dpng')