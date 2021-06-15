clc
clear all

addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/CONTACTMODELS/')
addpath('../ROUTINES/QUASISTATIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/QUADRATURE/')
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
s = 0.85;  H = 1.0;

%% Setup PCE in mu alone
Nq_pce = 10;
[xi, wi] = LAGWT(Nq_pce);
% wi = wi*(s/H);  % multiplied weight by mean
mui = xi*(s/H);
%%
tic
for ip=9:Nq_pce
    mu = mui(ip)*ones(MESH.Ne*MESH.Nq^2,1);

    %% Contact Model
    Estar = E/(1-nu^2)/2;
    Gstar = E/((1+nu)*(2-nu))/2;

    cn = Estar*sqrt(pi*Rad./lam.^3).*Nasps./Aels;
    ct = 4*Gstar*sqrt(pi*Rad./lam).*Nasps./Aels;

    lamo = lam;  muo = mu; gapo = gap;
    cno = cn.*exp(-lam.*z0); cto = ct.*exp(-lam.*z0);

    %% Create Object
    GM = MDOFGEN(M, K, zeros(size(M)), L);

    GM = GM.SETNLFUN(2+5, ...
        kron(MESH.Qm, eye(3))*L(1:(MESH.Nn*3), :), ...
        @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lamo, muo, gapo, varargin{:}), ...
        L(1:(MESH.Nn*3), :)'*kron(MESH.Tm, eye(3)));

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
%     load('rqnmsolmm1.mat', 'Ustat');
%     U0 = Ustat;
%     [Ustat, ~, eflag, ~, Jstat] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
%     if eflag<=0  % No convergence
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
%     end

    %% Linearized Analysis
    [Vstat, Wstat] = eigs(Jstat, M, 10, 'SM');
    [Wstat, si] = sort(sqrt(diag(Wstat)));
    Vstat = Vstat(:, si);
    Vstat = Vstat./sqrt(diag(Vstat'*M*Vstat))';

    %% March
    Na = 10;
    As = logspace(-0.5, 2.5, Na)*9.81/Wstat(1)^2;
    % As = logspace(-6, -3, Na);
    As = [-As(end:-1:1) As]';
    Eflags = zeros(1, 2*Na);
    load('rqnmsolmm1.mat', 'UlC')
    % UlC = zeros(GM.Ndofs+1, 2*Na);
    dUdalC = zeros(GM.Ndofs+1, 2*Na);

    ul0 = [Ustat+Vstat(:,1)*As(Na+1); Wstat(1)^2];
    fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
    fopts.Display = 'off';
    for ia=1:Na    
%         [UlC(:, Na+ia), ~, Eflags(Na+ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(As(Na+ia))], 1, Fv*Prestress, Ustat), UlC(:, Na+ia), fopts);
        [UlC(:, Na+ia), ~, Eflags(Na+ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(As(Na+ia))], 1, Fv*Prestress, Ustat), ul0, fopts);
        if eflag<=0
            [UlC(:, Na+ia), ~, Eflags(Na+ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(As(Na+ia))], 1, Fv*Prestress, Ustat), ul0, fopts);
            if eflag<=0
                error('No Convergence')
            end
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

    ul0 = [Ustat+Vstat(:,1)*As(Na+1-1); Wstat(1)^2];
    for ia=1:Na
%         [UlC(:, Na+1-ia), ~, Eflags(Na+1-ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(-As(Na+1-ia))], 1, Fv*Prestress, Ustat), UlC(:, Na+1-ia), fopts); 
        [UlC(:, Na+1-ia), ~, Eflags(Na+1-ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(-As(Na+1-ia))], 1, Fv*Prestress, Ustat), ul0, fopts);
        if eflag<=0
            [UlC(:, Na+1-ia), ~, Eflags(Na+1-ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(-As(Na+1-ia))], 1, Fv*Prestress, Ustat), ul0, fopts);
            if eflag<=0
                error('No Convergence')
            end
        end

        dRdA = [zeros(GM.Ndofs,1);-2*UlC(end, Na+1-ia)*As(Na+1-ia)];
        dUdalC(:, Na+1-ia) = -dRdUl\dRdA;
        if ia<Na
            ul0 = UlC(:, Na+1-ia);
            ul0(1:end-1) = Ustat+(ul0(1:end-1)-Ustat)*As(Na+1-ia-1)/As(Na+1-ia);
        end
        fprintf('%d, ', ia);
    end
    fprintf('\n');

    %% Post Processing (Hermite Interpolation + SHB)
    Ln = reshape([UlC(end,:); dUdalC(end,:)], 2*size(UlC,2), 1);

    Nq = 100;
    Qs = logspace(-5.5, -2.625, Nq)';

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
    Zts = zeros(Nq, 1);
    parfor (qi=1:Nq,8)
    % for qi=1:Nq
        Zts(qi) = GETFOURIERCOEFF(0, sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M +...
            squeeze(Udot(:, qi, :))*GM.C +...
            squeeze(Ut(:, qi, :))*GM.K + ...
            GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)), tol, Nits)),2))/(Qs(qi)^2*Lams(qi)^1.5);
        fprintf('%d\n', qi)
    end

    %% Save Information into file
    save(sprintf('./MUPCE/mupce_%d_N%d.mat', ip, Nq_pce), 'Qs', 'Phi', 'Lams', 'Zts', 'Wstat', 'Ustat');

    fprintf('=============================================\n')
    fprintf('Done %d/%d\n', ip, Nq_pce);
    fprintf('=============================================\n')
end
toc

%% Obtain PCE Coefficients
Nq_pce = 10;
[xi, wi] = LAGWT(Nq_pce);
Npce = Nq_pce-1;

[Psi, Integs] = PLAGU(0:Npce, xi);  %% Polynomials and their squared L2 norms
load(sprintf('./MUPCE/mupce_%d_N%d.mat', 1, 10), 'Qs')

Rcofs = zeros(length(Qs), Npce+1);
Wcofs = zeros(length(Qs), Npce+1);
Zcofs = zeros(length(Qs), Npce+1);
Wstatcofs = zeros(10, Npce+1);
for iq=1:Nq_pce
    load(sprintf('./MUPCE/mupce_%d_N%d.mat', iq, Nq_pce), 'Qs', 'Lams', 'Zts', 'Phi', 'Wstat')
    
    Rx = real((Qs.*(R(3,:)*Phi)').*Lams);
    Rcofs = Rcofs + wi(iq)*(Rx.*Psi(iq,:))./Integs';
    Wcofs = Wcofs + wi(iq)*(sqrt(Lams(:)).*Psi(iq, :))./Integs';
    Zcofs = Zcofs + wi(iq)*(Zts(:).*Psi(iq, :))./Integs';
    
    Wstatcofs = Wstatcofs + wi(iq)*(Wstat(:).*Psi(iq, :))./Integs';
end

%% Get PCE as polynomial coefficients with symbolics
x = sym('x');
assume(x, 'real')
rx = Rcofs*transpose(PLAGU(0:Npce, x));
wx = Wcofs*transpose(PLAGU(0:Npce, x));
zx = Zcofs*transpose(PLAGU(0:Npce, x));
wsx = Wstatcofs*transpose(PLAGU(0:Npce, x));
rxps = zeros(length(Qs), Npce+1);
wxps = zeros(length(Qs), Npce+1);
zxps = zeros(length(Qs), Npce+1);
wsxps = zeros(10, Npce+1);
for iq=1:length(Qs)
    rxps(iq, :) = sym2poly(rx(iq));
    wxps(iq, :) = sym2poly(wx(iq));
    zxps(iq, :) = sym2poly(zx(iq));
end
.Tick

save(sprintf('./MUPCE/mupce_N%d_cofs.mat', Nq_pce), 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'Integs', 'wxps', 'zxps', 'rxps', 'wsxps')

%% Confidence Intervals Over all points
W1s = zeros(size(Wcofs(:,1)));
W2s = zeros(size(Wcofs(:,1)));
Z1s = zeros(size(Zcofs(:,1)));
Z2s = zeros(size(Zcofs(:,1)));

alpha = 0.05;
opt = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
for iq=1:length(W1s)
    w10 = polyval(wxps(iq,:), -log(1-alpha/2));
    dwdx10 = polyval(polyder(wxps(iq,:)), -log(1-alpha/2));
    if sign(dwdx10)==-1
        w20 = w10;
        w10 = polyval(wxps(iq,:), -log(alpha/2));
    else
        w20 = polyval(wxps(iq,:), -log(alpha/2));
    end    
%     [W1s(iq), ~, eflag] = fsolve(@(w) POLYCDF(wxps(iq,:), w, ...
%                                               @(z) 1-exp(-z), ...
%                                               @(z) exp(-z), [0 inf], ...
%                                               alpha/2), w10, opt);
%     [W2s(iq), ~, eflag] = fsolve(@(w) POLYCDF(wxps(iq,:), w, ...
%                                               @(z) 1-exp(-z), ...
%                                               @(z) exp(-z), [0 inf], ...
%                                               1-alpha/2), w20, opt);
    [W1s(iq), ~, eflag] = fzero(@(w) POLYCDF(wxps(iq,:), w, ...
                                              @(z) 1-exp(-z), ...
                                              @(z) exp(-z), [0 inf], ...
                                              alpha/2), w10);
    [W2s(iq), ~, eflag] = fzero(@(w) POLYCDF(wxps(iq,:), w, ...
                                              @(z) 1-exp(-z), ...
                                              @(z) exp(-z), [0 inf], ...
                                              1-alpha/2), w20);

    z10 = polyval(zxps(iq,:), -log(1-alpha/2));
    dzdx10 = polyval(polyder(zxps(iq,:)), -log(1-alpha/2));
    if sign(dzdx10)==-1
        z20 = z10;
        z10 = polyval(zxps(iq,:), -log(alpha/2));
    else
        z20 = polyval(zxps(iq,:), -log(alpha/2));
    end
%     [Z1s(iq), ~, eflag] = fsolve(@(z) POLYCDF(zxps(iq,:), z, ...
%                                               @(z) 1-exp(-z), ...
%                                               @(z) exp(-z), [0 inf], ...
%                                               alpha/2), z10, opt); 
%     [Z2s(iq), ~, eflag] = fsolve(@(z) POLYCDF(zxps(iq,:), z, ...
%                                               @(z) 1-exp(-z), ...
%                                               @(z) exp(-z), [0 inf], ...
%                                               1-alpha/2), z20, opt);
    [Z1s(iq), ~, eflag] = fzero(@(z) POLYCDF(zxps(iq,:), z, ...
                                              @(z) 1-exp(-z), ...
                                              @(z) exp(-z), [0 inf], ...
                                              alpha/2), z10); 
    [Z2s(iq), ~, eflag] = fzero(@(z) POLYCDF(zxps(iq,:), z, ...
                                              @(z) 1-exp(-z), ...
                                              @(z) exp(-z), [0 inf], ...
                                              1-alpha/2), z20);
    
    fprintf('%d\n', iq);
end

%% 
load('./MATFILES/SaveFile_2021-Feb-19_High.mat', 'AMP_avg', 'FRE_avg', 'DAM_avg');
amin = 2.49;

figi = 2;

if any(arrayfun(@(x) x.Number, findall(0, 'type', 'figure'))==1)
    figure(figi); clf();
else
    figure(figi); clf();
    set(gcf, 'Position', [2800 550 1200 480])
end
set(gcf, 'Color', 'white')
for iq=1:Nq_pce
    load(sprintf('./MUPCE/mupce_%d_N%d.mat', iq, Nq_pce), 'Qs', 'Lams', 'Zts', 'Phi')

    figure(figi)
    subplot(1,2,1)
    ad=semilogx(Rx/9.81, sqrt(Lams)/2/pi, '-', 'Color', 0.6*[1 1 1]); hold on

    subplot(1,2,2)
    semilogx(Rx/9.81, Zts*100, '-', 'Color', 0.6*[1 1 1]); hold on
end
% Variances
Wvar = (Wcofs(:, 2:end).^2)*Integs(2:end);
Zvar = (Zcofs(:, 2:end).^2)*Integs(2:end);
Rvar = (Rcofs(:, 2:end).^2)*Integs(2:end);

subplot(1,2,1)
am = semilogx(Rcofs(:,1)/9.81, Wcofs(:,1)/2/pi, 'b-', 'LineWidth', 2);
% sm = fill(Rcofs([1:end end:-1:1],1)/9.81,
% [Wcofs(:,1)+3*sqrt(Wvar);
% Wcofs(end:-1:1,1)-3*sqrt(Wvar(end:-1:1))]/2/pi, 0.8*[1 1 1],
% 'EdgeColor', 'k', 'FaceAlpha', 0.4); hold on
sm = fill(Rcofs([1:end end:-1:1],1)/9.81, [W1s; W2s(end:-1:1)]/2/pi, ...
          0.8*[1 1 1], 'EdgeColor', 'k', 'FaceAlpha', 0.4); hold on
em = semilogx(AMP_avg(AMP_avg>=amin), FRE_avg(AMP_avg>=amin), 'k-', 'LineWidth', 2);
ll=legend([ad am sm em], 'Model Evaluations', 'Mean of PCE', ...
          sprintf('$%d \\%% $ C.I.', (1-alpha)*100), ...
          'Experimental Measurement', 'Location', 'southwest');
xlabel('Response Amplitude (g)')
ylabel('Natural Frequency (Hz)')
xlim([min(Rcofs(:,1)) max(Rcofs(:,1))]/9.81)

subplot(1,2,2)
semilogx(Rcofs(:,1)/9.81, Zcofs(:,1)*100, 'b-', 'LineWidth', 2)
% fill(Rcofs([1:end end:-1:1],1)/9.81, [Zcofs(:,1)+3*sqrt(Zvar);
% Zcofs(end:-1:1,1)-3*sqrt(Zvar(end:-1:1))]*100, 0.8*[1 1 1],
% 'EdgeColor', 'k', 'FaceAlpha', 0.4); hold on
fill(Rcofs([1:end end:-1:1],1)/9.81, [Z1s; Z2s(end:-1:1)]*100, ...
     0.8*[1 1 1], 'EdgeColor', 'k', 'FaceAlpha', 0.4); hold on
semilogx(AMP_avg(AMP_avg>=amin), DAM_avg(AMP_avg>=amin)*100, 'k-', 'LineWidth', 2);
xlabel('Response Amplitude (g)')
ylabel('Damping Factor (\%)')
xlim([min(Rcofs(:,1)) max(Rcofs(:,1))]/9.81)

sgtitle('PCE Results for Coefficient of Friction Variation')
% export_fig('./MUPCE/BBFIG.png', '-dpng')

% %% Single line CDF
% Nwss = 101;
% FYw = zeros(Nwss, 1);
% fYw = zeros(Nwss, 1);
% FYz = zeros(Nwss, 1);
% fYz = zeros(Nwss, 1);

% iq = 1;

% Ws = linspace(175, 195, Nwss)*2*pi;
% Zs = linspace(-1.6, -1.2, Nwss)*1e-6;
% for is=1:Nwss
%     [FYw(is), fYw(is)] = POLYCDF(wxps(iq,:), Ws(is), @(z) 1-exp(-z), @(z) exp(-z), [0 inf]);
%     [FYz(is), fYz(is)] = POLYCDF(zxps(iq,:), Zs(is), @(z) 1-exp(-z), @(z) exp(-z), [0 inf]);
%     fprintf('%d: %f %f\n', is, FYw(is), FYz(is));
% end

% figure(2); clf(); plot(Zs, FYz, '.-')

% %% PDF & CDF (Pg. 115 in Kobayashi, Mark & Turin)
% ps = linspace(0, 1, Nps+1); ps(end)=[];
% mxs = double(solve(diff(wx(iw))==0, x));
% 
% Wexts = double(subs(wx(iw), {x}, {mxs}));
% % mxs = mxs(Wexts>0);
% % Wexts = Wexts(Wexts>0);
% mxs = mxs(Wexts>0 & Wexts<220*2*pi);
% Wexts = Wexts(Wexts>0 & Wexts<220*2*pi);
% % Wminmax = minmax(Wexts(:)');
% Nwss = 100;
% Ws = linspace(min(Wexts), max(Wexts), Nwss);
% % fY = zeros(1, Nwss);
% FY = zeros(1, Nwss);
% 
% parfor is=1:Nwss
%     x_rts = double(solve(wx(iw)==Ws(is), x));
% %     if sum(imag(x_rts)==0 & x_rts>=0)==0
% %         keyboard
% %     end
%     x_rts = x_rts(imag(x_rts)==0 & x_rts>=0);  % only real roots > 0
%     der_x = double(subs(diff(wx(iw)), {x}, {x_rts}));
%     fY(is) = sum(exp(-x_rts)./abs(der_x));
%     if length(x_rts)==1
%         FY(is) = 1-exp(-x_rts);
%     else
%         FY(is) = sum((1-exp(-x_rts)).*sign(der_x));
%     end
%     fprintf('%d\n', is);
% end

%     %%
% tic
% fsolve(@(w) POLYCDF(sym2poly(wx(iw)), w, @(z) 1-exp(-z), @(z) exp(-z), alpha/2), 170*2*pi, opt)
% toc
% %%
% figure(2)
% clf()
% plot(Ws/2/pi, FY, '.-'); hold on
% % plot(Ws/2/pi, FY,'.-')
% plot(xlim, alpha/2*[1 1], 'k-')
% plot(xlim, (1-alpha/2)*[1 1], 'k-')
% plot(W1*[1 1]/2/pi, ylim, 'k-')
% plot(W2*[1 1]/2/pi, ylim, 'k-')
% %%
% figure(3)
% clf()
% % FYc = cumtrapz(Ws, fY); FYc(FYc>1) = 1;
% % plot(Ws, FYc); hold on; 
% plot(Ws, FY); hold on
% Nsamp = 1001;
% rng(1)
% Wsamps = Wcofs*PLAGU(0:Npce, exprnd(1, Nsamp,1))';
% [ps, ws] = ecdf(Wsamps(iw,:)); 
% plot(ws, ps, '-'); grid on
% 
% % FYs = zeros(size(ws));
% % for iw=1:length(ws)
% %     FYs(iw) = Fcdf(ws(iw));
% % end
