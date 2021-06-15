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

% Collect
Nasps = kron(Nasps, ones(MESH.Nq^2,1));
lam   = kron(lam, ones(MESH.Nq^2,1));
z0    = log(Nasps)./lam;

%% Curvature Radii
R1 = (R1top.CRAD(:,1).*R1top.NASPS(:,1)+R1bot.CRAD(:,1).*R1bot.NASPS(:,1))./(R1top.NASPS(:,1)+R1bot.NASPS(:,2));
R2 = (R2top.CRAD(:,1).*R2top.NASPS(:,1)+R2bot.CRAD(:,1).*R2bot.NASPS(:,1))./(R2top.NASPS(:,1)+R2bot.NASPS(:,2));

Rad = (R1+R2)/2;
Rad = kron(Rad, ones(MESH.Nq^2,1));

%% Prestress
Prestress = mean([12002 12075 12670]);

%% Coefficient of Friction (AISI 304N SS)
% s = 186e6; H = 655e6;
% s = 190e6;  H = 545e6;
s = 0.85;  H = 1.0;

%% Setup PCE in microscale params alone
Nq_pce = 10;
[xi, wi] = GPHWT(Nq_pce);

lamt1i = R1top.LLX0s_sd(:,1)+R1top.LLX0s_sd(:,3).*xi(:)';
lamt2i = R2top.LLX0s_sd(:,1)+R2top.LLX0s_sd(:,3).*xi(:)';
lamb1i = R1bot.LLX0s_sd(:,1)+R1bot.LLX0s_sd(:,3).*xi(:)';
lamb2i = R2bot.LLX0s_sd(:,1)+R2bot.LLX0s_sd(:,3).*xi(:)';
lam1i = (R1top.NASPS(:,1)+R1bot.NASPS(:,1))./(R1top.NASPS(:,1)./lamt1i+R1bot.NASPS(:,1)./lamb1i);
lam2i = (R2top.NASPS(:,1)+R2bot.NASPS(:,1))./(R2top.NASPS(:,1)./lamt2i+R2bot.NASPS(:,1)./lamb2i);
lami = (lam1i+lam2i)/2;

lami = kron(lami, ones(Nq^2,1));
z0i = log(Nasps)./lami;
clear lamt1i lamt2i lamb1i lamb2i lam1i lam2i

%%
tic
for ip=1:Nq_pce
    mu = (s/H)*ones(MESH.Ne*MESH.Nq^2,1);

    %% Contact Model
    Estar = E/(1-nu^2)/2;
    Gstar = E/((1+nu)*(2-nu))/2;

    lamo = lami(:,ip); muo = mu; gapo = gap; z0o = z0i(:,ip);

    cn = Estar*sqrt(pi*Rad./lamo.^3).*Nasps./Aels;
    ct = 4*Gstar*sqrt(pi*Rad./lamo).*Nasps./Aels;

    cno = cn.*exp(-lamo.*z0o); cto = ct.*exp(-lamo.*z0o);

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
    As = [-As(end:-1:1) As]';
    Eflags = zeros(1, 2*Na);
    UlC = zeros(GM.Ndofs+1, 2*Na);
    dUdalC = zeros(GM.Ndofs+1, 2*Na);

    ul0 = [Ustat+Vstat(:,1)*As(Na+1); Wstat(1)^2];
    fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
    fopts.Display = 'off';
    for ia=1:Na    
        [UlC(:, Na+ia), ~, Eflags(Na+ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(As(Na+ia))], 1, Fv*Prestress, Ustat), ul0, fopts);
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

    ul0 = [Ustat+Vstat(:,1)*As(Na+1-1); Wstat(1)^2];
    for ia=1:Na
        [UlC(:, Na+1-ia), ~, Eflags(Na+1-ia), ~, dRdUl] = fsolve(@(ul) GM.RQMRESFUN([ul; log10(-As(Na+1-ia))], 1, Fv*Prestress, Ustat), ul0, fopts);
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
    save(sprintf('./MSCPCE/mscpce_%d_N%d.mat', ip, Nq_pce), 'Qs', 'Phi', 'Lams', 'Zts', 'Wstat', 'Ustat');

    fprintf('=============================================\n')
    fprintf('Done %d/%d\n', ip, Nq_pce);
    fprintf('=============================================\n')
end
toc

%% Obtain PCE Coefficients
Nq_pce = 10;
[xi, wi] = GPHWT(Nq_pce);
Npce = Nq_pce-1;

[Psi, Integs] = PHERM(0:Npce, xi);  %% Polynomials and their squared L2 norms
load(sprintf('./MSCPCE/mscpce_%d_N%d.mat', 1, 10), 'Qs')

Rcofs = zeros(length(Qs), Npce+1);
Wcofs = zeros(length(Qs), Npce+1);
Zcofs = zeros(length(Qs), Npce+1);
Wstatcofs = zeros(10, Npce+1);
for iq=1:Nq_pce
    load(sprintf('./MSCPCE/mscpce_%d_N%d.mat', iq, Nq_pce), 'Qs', 'Lams', 'Zts', 'Phi', 'Wstat')
    
    Rx = real((Qs.*(R(3,:)*Phi)').*Lams);
    Rcofs = Rcofs + wi(iq)*(Rx.*Psi(iq,:))./Integs';
    Wcofs = Wcofs + wi(iq)*(sqrt(Lams(:)).*Psi(iq, :))./Integs';
    Zcofs = Zcofs + wi(iq)*(Zts(:).*Psi(iq, :))./Integs';
    
    Wstatcofs = Wstatcofs + wi(iq)*(Wstat(:).*Psi(iq, :))./Integs';
end

%% Get PCE as polynomial coefficients with symbolics
x = sym('x');
assume(x, 'real')
rx = Rcofs*transpose(PHERM(0:Npce, x));
wx = Wcofs*transpose(PHERM(0:Npce, x));
zx = Zcofs*transpose(PHERM(0:Npce, x));
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
for iw=1:10
    wsxps(iw, :) = sym2poly(wsx(iw));
end


save(sprintf('./MSCPCE/mscpce_N%d_cofs.mat', Nq_pce), 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'Integs', 'wxps', 'zxps', 'rxps', 'wsxps')

%% Confidence Intervals Over all points
W1s = zeros(size(Wcofs(:,1)));
W2s = zeros(size(Wcofs(:,1)));
Z1s = zeros(size(Zcofs(:,1)));
Z2s = zeros(size(Zcofs(:,1)));

alpha = 0.05;
opt = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
for iq=1:length(W1s)
    w10 = polyval(wxps(iq,:), norminv(alpha/2));
    dwdx10 = polyval(polyder(wxps(iq,:)), norminv(alpha/2));
    if sign(dwdx10)==-1
        w20 = w10;
        w10 = polyval(wxps(iq,:), norminv(1-alpha/2));
    else
        w20 = polyval(wxps(iq,:), norminv(1-alpha/2));
    end
    [W1s(iq), ~, eflag] = fsolve(@(w) POLYCDF(wxps(iq,:), w, ...
                                              @(z) normcdf(z), ...
                                              @(z) normpdf(z), [-inf inf], ...
                                              alpha/2), w10, opt);
    [W2s(iq), ~, eflag] = fsolve(@(w) POLYCDF(wxps(iq,:), w, ...
                                              @(z) normcdf(z), ...
                                              @(z) normpdf(z), [-inf inf], ...
                                              1-alpha/2), w20, opt);

    z10 = polyval(zxps(iq,:), norminv(alpha/2));
    dzdx10 = polyval(polyder(zxps(iq,:)), norminv(alpha/2));
    if sign(dzdx10)==-1
        z20 = z10;
        z10 = polyval(zxps(iq,:), norminv(1-alpha/2));
    else
        z20 = polyval(zxps(iq,:), norminv(1-alpha/2));
    end
    [Z1s(iq), ~, eflag] = fsolve(@(z) POLYCDF(zxps(iq,:), z, ...
                                              @(z) normcdf(z), ...
                                              @(z) normpdf(z), [-inf inf], ...
                                              alpha/2), z10, opt); 
    [Z2s(iq), ~, eflag] = fsolve(@(z) POLYCDF(zxps(iq,:), z, ...
                                              @(z) normcdf(z), ...
                                              @(z) normpdf(z), [-inf inf], ...
                                              1-alpha/2), z20, opt);
    
    fprintf('%d\n', iq);
end

%% Process PCE & make Plots
load('./MATFILES/SaveFile_2021-Feb-19_High.mat', 'AMP_avg', 'FRE_avg', 'DAM_avg');
amin = 2.49;

if any(arrayfun(@(x) x.Number, findall(0, 'type', 'figure'))==1)
    figure(1); clf();
else
    figure(1); clf();
    set(gcf, 'Position', [2800 550 1200 480])
end
set(gcf, 'Color', 'white')
for iq=1:Nq_pce
    load(sprintf('./MSCPCE/mscpce_%d_N%d.mat', iq, Nq_pce), 'Qs', 'Lams', 'Zts', 'Phi')

    figure(1)
    subplot(1,2,1)
    ad=semilogx(Rx/9.81, sqrt(Lams)/2/pi, '-', 'Color', 0.6*[1 1 1]); hold on

    subplot(1,2,2)
    semilogx(Rx/9.81, Zts*100, '-', 'Color', 0.6*[1 1 1]); hold on
end
Wvar = (Wcofs(:, 2:end).^2)*Integs(2:end);
Zvar = (Zcofs(:, 2:end).^2)*Integs(2:end);
Rvar = (Rcofs(:, 2:end).^2)*Integs(2:end);

save(sprintf('./MSCPCE/mscpce_N%d_cofs.mat', Nq_pce), 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Integs')

subplot(1,2,1)
am = semilogx(Rcofs(:,1)/9.81, Wcofs(:,1)/2/pi, 'b-', 'LineWidth', 2);
sm = fill(Rcofs([1:end end:-1:1],1)/9.81, [Wcofs(:,1)+3*sqrt(Wvar); Wcofs(end:-1:1,1)-3*sqrt(Wvar(end:-1:1))]/2/pi, 0.8*[1 1 1], 'EdgeColor', 'k', 'FaceAlpha', 0.4); hold on
% sm = fill(Rcofs([1:end end:-1:1],1)/9.81, [Wcofs(:,1)+3*sqrt(Wvar); Wcofs(end:-1:1,1)-3*sqrt(Wvar(end:-1:1))]/2/pi, 0.8*[1 1 1], 'EdgeColor', 'k', 'FaceAlpha', 0.4); hold on
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

sgtitle('PCE Results for Asperity Mean Height Variation')
export_fig('./MSCPCE/BBFIG.png', '-dpng')