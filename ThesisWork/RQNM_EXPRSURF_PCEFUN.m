function [] = RQNM_EXPRSURF_PCEFUN(Ixs, nxi, Nq_pces, pref)
%RQNM_EXPRSURF_PCEFUN conducts PCE evaluation for given Design
% Order of parameters Assumed: [mu, msc, prestress, rotx, roty, gap]

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
    Nasps1 = R1top.NASPS(:,1)+R1bot.NASPS(:,1);
    Nasps2 = R2top.NASPS(:,1)+R2bot.NASPS(:,1);
    Nasps = (Nasps1+Nasps2)/2;

    lam1 = (R1top.NASPS(:,1)+R1bot.NASPS(:,1))./(R1top.NASPS(:,1)./R1top.LLX0s_sd(:,1)+R1bot.NASPS(:,1)./R1bot.LLX0s_sd(:,1));
    lam2 = (R2top.NASPS(:,1)+R2bot.NASPS(:,1))./(R2top.NASPS(:,1)./R2top.LLX0s_sd(:,1)+R2bot.NASPS(:,1)./R2bot.LLX0s_sd(:,1));
    lam = (lam1+lam2)/2;

    % Collect
    Nasps = kron(Nasps, ones(Nq^2,1));
    lam   = kron(lam, ones(Nq^2,1));
    % z0    = kron(z0, ones(Nq^2,1));
    % z0 = -log(0.1)./lam;
    z0 = log(Nasps)./lam;

    %% Curvature Radii
    R1 = (R1top.CRAD(:,1).*R1top.NASPS(:,1)+R1bot.CRAD(:,1).*R1bot.NASPS(:,1))./(R1top.NASPS(:,1)+R1bot.NASPS(:,1));
    R2 = (R2top.CRAD(:,1).*R2top.NASPS(:,1)+R2bot.CRAD(:,1).*R2bot.NASPS(:,1))./(R2top.NASPS(:,1)+R2bot.NASPS(:,1));

    Rad = (R1+R2)/2;
    Rad = kron(Rad, ones(Nq^2,1));

%     %% Setup PCE in mu alone
%     Nq_pce = 10;
%     [xi, wi] = LAGWT(Nq_pce);
%     % wi = wi*(s/H);  % multiplied weight by mean
%     mui = xi*(s/H);

    %% Use the Given PCE Coefficients: [mu, msc, prestress, rotx, roty, gap]
    Xis = zeros(6,1);
    Wis = zeros(6,1);
    
    % 1. Coefficient of Friction
    % s = 190e6;  H = 545e6; (AISI 304N SS)
    s = 0.85;  H = 1.0;
    [xi, wi] = LAGWT(Nq_pces(1));
    mu = xi(Ixs(1))*(s/H)*ones(MESH.Ne*MESH.Nq^2, 1);
    Xis(1) = xi(Ixs(1));  Wis(1) = wi(Ixs(1));
    
    % 2. Micro-Scale Lambda
    [xi, wi] = GPHWT(Nq_pces(2));
    lamt1i = R1top.LLX0s_sd(:,1)+R1top.LLX0s_sd(:,3)*xi(Ixs(2));
    lamt2i = R2top.LLX0s_sd(:,1)+R2top.LLX0s_sd(:,3)*xi(Ixs(2));
    lamb1i = R1bot.LLX0s_sd(:,1)+R1bot.LLX0s_sd(:,3)*xi(Ixs(2));
    lamb2i = R2bot.LLX0s_sd(:,1)+R2bot.LLX0s_sd(:,3)*xi(Ixs(2));
    lam1i = (R1top.NASPS(:,1)+R1bot.NASPS(:,1))./(R1top.NASPS(:,1)./lamt1i+R1bot.NASPS(:,1)./lamb1i);
    lam2i = (R2top.NASPS(:,1)+R2bot.NASPS(:,1))./(R2top.NASPS(:,1)./lamt2i+R2bot.NASPS(:,1)./lamb2i);
    lam = (lam1i+lam2i)/2;
    
    lam = kron(lam, ones(Nq^2, 1));
    z0 = log(Nasps)./lam;
    clear lamt1i lamt2i lamb1i lamb2i lam1i lam2i
    Xis(2) = xi(Ixs(2));  Wis(2) = wi(Ixs(2));

    % 3. Prestress
    Plevels = [12002 12075 12670];
    Psd = 2100;
%     Psd = 600;
    [xi, wi] = GPHWT(Nq_pces(3));
    Prestress = mean(Plevels)+Psd*xi(Ixs(3));
    Xis(3) = xi(Ixs(3));  Wis(3) = wi(Ixs(3));
    
    % 4-5. Rotx-Roty
    [xi1, wi1] = GPHWT(Nq_pces(4));
    [xi2, wi2] = GPHWT(Nq_pces(5));
    theta_sd = deg2rad(15);
    thetas = theta_sd*[xi1(Ixs(4)); xi2(Ixs(5))];
    % 6. Gap
    [xi, wi] = GPHWT(Nq_pces(6));

    TFM = [1,  0               , 0; 
           0,  cos(thetas(1)), sin(thetas(1));
           0, -sin(thetas(1)), cos(thetas(1))]*...
        [cos(thetas(2)), 0, -sin(thetas(2));
         0               , 1,  0;
         sin(thetas(2)), 0,  cos(thetas(2))];  % Transformation for gap
     gtops = [(R1top.BilinPlaneQPs(:,1)+xi(Ixs(6))*R1top.BilinPlaneQPs(:,2)) (R2top.BilinPlaneQPs(:,1)+xi(Ixs(6))*R2top.BilinPlaneQPs(:,2))];
     gbots = [(R1bot.BilinPlaneQPs(:,1)+xi(Ixs(6))*R1bot.BilinPlaneQPs(:,2)) (R2bot.BilinPlaneQPs(:,1)+xi(Ixs(6))*R2bot.BilinPlaneQPs(:,2))];
     xygs = [MESH.Qm*MESH.Nds gtops gbots];
     gap1r = xygs(:, [1 2 3])*TFM(:, 3) - xygs(:, [1 2 5])*TFM(:, 3);
     gap2r = xygs(:, [1 2 4])*TFM(:, 3) - xygs(:, [1 2 6])*TFM(:, 3);
     
     gapr = (gap1r+gap2r)/2;
     
     Xis(4:5) = [xi1(Ixs(4));xi2(Ixs(5))];  Wis(4:5) = [wi1(Ixs(4));wi2(Ixs(5))];
     Xis(6) = xi(Ixs(6));  Wis(6) = wi(Ixs(6));
     
    %% Run Simulations
    tic
    %% Contact Model
    Estar = E/(1-nu^2)/2;
    Gstar = E/((1+nu)*(2-nu))/2;

    cn = Estar*sqrt(pi*Rad./lam.^3).*Nasps./Aels;
    ct = 4*Gstar*sqrt(pi*Rad./lam).*Nasps./Aels;

    cno = cn.*exp(-lam.*z0); cto = ct.*exp(-lam.*z0);
    %% Create Object
    GM = MDOFGEN(M, K, zeros(size(M)), L);

    GM = GM.SETNLFUN(2+5, ...
        kron(MESH.Qm, eye(3))*L(1:(MESH.Nn*3), :), ...
        @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lam, mu, gap, varargin{:}), ...
        L(1:(MESH.Nn*3), :)'*kron(MESH.Tm, eye(3)));

    %% Linearized Stiffness
    knlin = lam.*Prestress/Aint;
    ktlin = cto./cno*Prestress/Aint;

    K0 = zeros(size(L,1));
    K0(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = MESH.Tm*diag(ktlin)*MESH.Qm;
    K0(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = MESH.Tm*diag(ktlin)*MESH.Qm;
    K0(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = MESH.Tm*diag(knlin)*MESH.Qm;
    K0 = L'*K0*L;

    %% Prestress Analysis
    fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');

    U0 = (K+K0)\(Fv*Prestress);
    U0 = L*U0;  U0(3:3:MESH.Nn*3) = MESH.Qm\gap;
    U0 = L\U0;

    GM.NLTs.func = @(t, u, varargin) EXPROUGHFRICT(t, u, cto*0, cno, lam, mu, gap, varargin{:});
    [U0, ~, eflag] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
    if eflag<=0
        load('rqnmsolmm1.mat', 'Ustat');
        U0 = Ustat;
        U0 = L*U0;  U0(3:3:MESH.Nn*3) = MESH.Qm\gap;
        U0 = L\U0;

        [U0, ~, eflag] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);
        if eflag<=0
            error('No Prestress Convergence')
        end
    end
    GM.NLTs.func = @(t, u, varargin) EXPROUGHFRICT(t, u, cto, cno, lam, mu, gap, varargin{:});
    [Ustat, ~, eflag, ~, Jstat] = fsolve(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, fopts);

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
    % load('rqnmsolmm1.mat', 'UlC')
    UlC = zeros(GM.Ndofs+1, 2*Na);
    dUdalC = zeros(GM.Ndofs+1, 2*Na);

    ul0 = [Ustat+Vstat(:,1)*As(Na+1); Wstat(1)^2];
    fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
    fopts.Display = 'off';
    for ia=1:Na
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
    save(sprintf('./ALLPCE/%s_%d.mat', pref, nxi), 'Qs', 'Phi', 'Lams', 'Zts', 'Wstat', 'Ustat', 'Xis', 'Wis');

    fprintf('=============================================\n')
    fprintf('Done %d\n', nxi);
    fprintf('=============================================\n')
    %%
    toc
end
