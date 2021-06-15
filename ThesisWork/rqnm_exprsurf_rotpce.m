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

%% Setup PCE in gap rotation alone
Nq_pce = 10;
[xi, wi] = GPHWT(Nq_pce);
theta_sd = deg2rad(15);
thetas = theta_sd*xi;

% xygs = [MESH.Qm*MESH.Nds gap];
xygs = [MESH.Qm*MESH.Nds R1top.BilinPlaneQPs(:,1) R2top.BilinPlaneQPs(:,1) R1bot.BilinPlaneQPs(:,1) R2bot.BilinPlaneQPs(:,1)];

%%
for np=1:Nq_pce^2
    ip = fix((np-1)/Nq_pce);
    jp = (np-1)-ip*Nq_pce;
    
    TFM = [1,  0               , 0; 
           0,  cos(thetas(ip+1)), sin(thetas(ip+1));
           0, -sin(thetas(ip+1)), cos(thetas(ip+1))]*...
        [cos(thetas(jp+1)), 0, -sin(thetas(jp+1));
         0               , 1,  0;
         sin(thetas(ip+1)), 0,  cos(thetas(ip+1))];  % Transformation for gap
    gap1r = xygs(:, [1 2 3])*TFM(:, 3) - xygs(:, [1 2 5])*TFM(:, 3);
    gap2r = xygs(:, [1 2 4])*TFM(:, 3) - xygs(:, [1 2 6])*TFM(:, 3);
    gapr = (gap1r+gap2r)/2;
    
    mu = (s/H)*ones(MESH.Ne*MESH.Nq^2,1);

    %% Contact Model
    Estar = E/(1-nu^2)/2;
    Gstar = E/((1+nu)*(2-nu))/2;

    cn = Estar*sqrt(pi*Rad./lam.^3).*Nasps./Aels;
    ct = 4*Gstar*sqrt(pi*Rad./lam).*Nasps./Aels;

    lamo = lam; muo = mu; gapo = gapr;
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
    As = [-As(end:-1:1) As]';
    Eflags = zeros(1, 2*Na);
    UlC = zeros(GM.Ndofs+1, 2*Na);
    dUdalC = zeros(GM.Ndofs+1, 2*Na);

    tic
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
    toc

    %% Post Processing (Hermite Interpolation + SHB)
    Ln = reshape([UlC(end,:); dUdalC(end,:)], 2*size(UlC,2), 1);

    Nq = 100;
    Qs = logspace(-5.5, -2.625, Nq)';

    Nt = 2^7;
    t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
    qt = cos(t).*Qs';
    tic
    [Lt, Nint, dNint] = HERMINTERP(As, Ln, qt(:));
    toc
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
        Zts(qi) = GETFOURIERCOEFF(0, sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M +...
            squeeze(Udot(:, qi, :))*GM.C +...
            squeeze(Ut(:, qi, :))*GM.K + ...
            GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)), tol, Nits)),2))/(Qs(qi)^2*Lams(qi)^1.5);
        fprintf('%d\n', qi)
    end
    toc

    %% Save Information into file
    save(sprintf('./ROTPCE/rotpce_%d_N%d.mat', np, Nq_pce), 'xi', 'wi', 'ip', 'jp', 'Qs', 'Phi', 'Lams', 'Zts', 'Wstat', 'Ustat');

    fprintf('=============================================\n')
    fprintf('Done %d/%d\n', np, Nq_pce^2);
    fprintf('=============================================\n')
end
return

%% Process PCE & make Plots
Nq_pce = 10;
[xi, wi] = GPHWT(Nq_pce);
theta_sd = deg2rad(15);
thetas = theta_sd*xi;

Npce = 9

IJs = [repmat((0:Npce)', Npce+1,1) kron((0:Npce)', ones(Npce+1,1))];
IJs = IJs(sum(IJs,2)<=Npce,:);
Npsi = size(IJs,1);

Psi = zeros(Nq_pce^2, Npsi);
Integs = zeros(Npsi, 1);
for np=1:Npsi
    [p1, i1] = PHERM(IJs(np,1), xi);
    [p2, i2] = PHERM(IJs(np,2), xi);
    
    Psi(:, np) = kron(p1, p2);
    Integs(np) = i1*i2;
end
load('./ROTPCE/rotpce_1_N10.mat', 'Qs')
Rcofs = zeros(length(Qs), Npsi);
Wcofs = zeros(length(Qs), Npsi);
Zcofs = zeros(length(Qs), Npsi);
Wstatcofs = zeros(10, Npsi);
Wss = zeros(Nq_pce, Nq_pce);
for nq=1:Nq_pce^2
    ip = fix((nq-1)/Nq_pce);
    jp = (nq-1)-ip*Nq_pce;

    load(sprintf('./ROTPCE/rotpce_%d_N%d.mat', nq, Nq_pce), 'Qs', 'Lams', 'Zts', 'Phi', 'Wstat')
    Wss(ip+1, jp+1) = Wstat(1);
    
    Rx = real((Qs.*(R(3,:)*Phi)').*Lams);
    Rcofs = Rcofs + wi(ip+1)*wi(jp+1)*(Rx.*Psi(nq,:))./Integs';
    Wcofs = Wcofs + wi(ip+1)*wi(jp+1)*(sqrt(Lams(:)).*Psi(nq, :))./Integs';
    Zcofs = Zcofs + wi(ip+1)*wi(jp+1)*(Zts(:).*Psi(nq, :))./Integs';
    Wstatcofs = Wstatcofs + wi(ip+1)*wi(jp+1)*(Wstat(:).*Psi(nq, :))./Integs';
end

save(sprintf('./ROTPCE/rotpce_N%d_cofs.mat', Nq_pce), 'Qs', 'Rcofs', 'Wcofs', 'Zcofs', 'Wstatcofs', 'Integs', 'IJs')
%% Plotting
load('./MATFILES/SaveFile_2021-Feb-19_High.mat', 'AMP_avg', 'FRE_avg', 'DAM_avg');
amin = 2.49;

if any(arrayfun(@(x) x.Number, findall(0, 'type', 'figure'))==1)
    figure(1); clf();
else
    figure(1); clf();
    set(gcf, 'Position', [2800 550 1200 480])
end
set(gcf, 'Color', 'white')
for nq=1:Nq_pce^2
    
    figure(1)
    subplot(1,2,1)
    ad=semilogx(Rx/9.81, sqrt(Lams)/2/pi, '-', 'Color', 0.6*[1 1 1]); hold on

    subplot(1,2,2)
    semilogx(Rx/9.81, Zts*100, '-', 'Color', 0.6*[1 1 1]); hold on
    
%     sgtitle(sprintf('($\\theta_x$, $\\theta_y$) = (%.2f, %.2f) degs', rad2deg(thetas(ip+1)), rad2deg(thetas(jp+1))));
end
Wvar = (Wcofs(:, 2:end).^2)*Integs(2:end);
Zvar = (Zcofs(:, 2:end).^2)*Integs(2:end);
Rvar = (Rcofs(:, 2:end).^2)*Integs(2:end);
Wstatvar = (Wstatcofs(:, 2:end).^2)*Integs(2:end);

subplot(1,2,1)
sm = fill(Rcofs([1:end end:-1:1],1)/9.81, [Wcofs(:,1)+sqrt(Wvar)*3; Wcofs(end:-1:1,1)-sqrt(Wvar(end:-1:1))*3]/2/pi, 0.8*[1 1 1], 'EdgeColor', 'k', 'FaceAlpha', 0.4); hold on
am = semilogx(Rcofs(:,1)/9.81, Wcofs(:,1)/2/pi, 'b-', 'LineWidth', 2);
em = semilogx(AMP_avg(AMP_avg>=amin), FRE_avg(AMP_avg>=amin), 'k-', 'LineWidth', 2);
legend([ad am sm em], 'Model Evaluations', 'Mean of PCE', 'PCE Mean $\pm$ 3 Std. Dev.', 'Experimental Measurement', 'Location', 'southwest')
xlabel('Response Amplitude (g)')
ylabel('Natural Frequency (Hz)')
xlim([min(Rcofs(:,1)) max(Rcofs(:,1))]/9.81)

subplot(1,2,2)
fill(Rcofs([1:end end:-1:1],1)/9.81, [Zcofs(:,1)+sqrt(Zvar)*3; Zcofs(end:-1:1,1)-sqrt(Zvar(end:-1:1))*3]*100, 0.8*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on
plot(Rcofs([1:end end:-1:1],1)/9.81, [Zcofs(:,1)+sqrt(Zvar)*3; Zcofs(end:-1:1,1)-sqrt(Zvar(end:-1:1))*3]*100, 'k-')
semilogx(Rcofs(:,1)/9.81, Zcofs(:,1)*100, 'b-', 'LineWidth', 2)
semilogx(AMP_avg(AMP_avg>=amin), DAM_avg(AMP_avg>=amin)*100, 'k-', 'LineWidth', 2);
xlabel('Response Amplitude (g)')
ylabel('Damping Factor (\%)')
xlim([min(Rcofs(:,1)) max(Rcofs(:,1))]/9.81)

sgtitle('PCE Results for Stage Attitude Variation')
% export_fig('./ROTPCE/BBFIG.png', '-dpng')