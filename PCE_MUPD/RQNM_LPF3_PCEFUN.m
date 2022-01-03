function [outputs] = RQNM_LPF3_PCEFUN(Ixs, nxi, parprops, Nq_pces, pref, varargin)
%RQNM_LPENCONF_PCEFUN conducts PCE evaluation for given design.
% Order of parameters assumed : [kts, kns, mus, gaps]
%
%   INPUTS :
%       Ixs
%       nxi
%       Nq_pces  : [nq_kts, nq_kns, nq_mus, nq_gaps]
% 	parprops : Parameter properties array of structures containing,
% 		L_pars   : Ls (Ne*Nq^2, num_random) for each parameter
% 		pdists 	 : makedist object
%    		quadfun  : function handle to generate quadr. pts & wts
%       pref
%%%%%%%% OPTIONAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       mdi         : Mode of interest (1, 2, 3)
%       range       : [AMIN, AMAX] Amplitude Range
%       simmode     : Simulation type (how to interpret 'Ixs')
%           'quad'      : As indices of quadrature points
%           'proper'    : As actual random number

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Choose Simulation Configurations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Choose Mode
    mds = [1 3 5];
    
    mdi = 1;
    if length(varargin)>=1
        mdi = varargin{1};
    end
    
    %% Choose Amplitude Range
    Arange = [-7.0, -3.0];
    if length(varargin)>=2
        Arange = varargin{2}(:)';
    end
    Qrange = Arange + [0.1 -0.1];
    
    %% Choose Simulation Mode
    simmode = 'quad';
    if length(varargin)>=3
        simmode = varargin{3};
    end
    
    %% Load Model
    model = 'BRB_Thesis';
    E = 1.9231e11;
    nu = 0.3;
    
    load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'K', 'M', ...
         'R', 'L', 'Fv', 'TFM');
    
    Nds = dlmread(sprintf('../MODELS/%s/Nodes.dat', model));
    Quad = dlmread(sprintf('../MODELS/%s/Elements.dat', model));

    Nq = 2;
    MESH = MESH2D(Nds, 3, [], Quad, Nq);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Prepare Contact Model Parameters (FIX THIS TO HAVE 1 value per q.p.) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    kt = zeros(1, MESH.Ne*MESH.Nq^2);
    kn = zeros(1, MESH.Ne*MESH.Nq^2);
    mu = zeros(1, MESH.Ne*MESH.Nq^2);
    gap = zeros(1, MESH.Ne*MESH.Nq^2);
%     parslist = {'kt', 'kn', 'mu', 'gap'};

    Xis = zeros(length(parprops), 1);
    Wis = zeros(length(parprops), 1);
    for pi=1:length(parprops)
        [xi, wi] = parprops(pi).quadfun(Nq_pces(pi));
        
        if strcmp(simmode, 'quad')
            Xis(pi) = xi(Ixs(pi));  Wis(pi) = wi(Ixs(pi));
            
            eval([parprops(pi).par '(:) = ' num2str(parprops(pi).map(Xis(pi)),'%.100e') ';']);
        else
            Xis(pi) = Ixs(pi);  Wis(pi) = 1.0;
            
            eval([parprops(pi).par '(:) = ' num2str(Xis(pi),'%.100e') ';']);
        end
    end
    Kt = [kt; kt; zeros(1, MESH.Ne*MESH.Nq^2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. Create MDOFGEN MODEL & Assign Contact Model %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GM = MDOFGEN(M, K, zeros(size(M)), L);
    
    fnl = @(t, u, varargin) ELDRYFRICT(t, u, Kt, kn, mu, gap, ...
                                       varargin{:});
    GM = GM.SETNLFUN(2+5, ...
                     kron(MESH.Qm, eye(3))*L(1:(MESH.Nn*3), :), ...
                     fnl, ...
                     L(1:(MESH.Nn*3), :)'*kron(MESH.Tm, eye(3)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % 4. START SIMULATIONS %
    %%%%%%%%%%%%%%%%%%%%%%%%
    tic
    %% Prestress Analysis
    opts = struct('reletol', 1e-6, 'Display', true, 'lsrch', 0);

    Plevels = [12002 12075 12670];
    Prestress = Plevels(1);
    
    % Initial Guess Fully Engaged Stiffness
    Kmat = [kt(1) 0 0; 0 kt(1) 0; 0 0 kn(1)];
    K0 = L(1:MESH.Nn*3,:)'*kron(MESH.Tm,eye(3))*...
         kron(eye(MESH.Ne*MESH.Nq^2),Kmat)*...
         kron(MESH.Qm, eye(3))*L(1:MESH.Nn*3,:);    
    U0 = (K+K0)\(Fv*Prestress);
    
    [Ustat, ~, eflag, ~, Jstat] = NSOLVE(@(U) GM.STATRESFUN(U, Fv*Prestress), ...
                                     U0, opts); 
    if eflag<=0
        load('statsolguess.mat', 'Ustat');
        U0 = Ustat;
        U0 = L*U0;  U0(3:3:MESH.Nn*3) = MESH.Qm\gap;
        U0 = L\U0;
        
        [Ustat, ~, eflag, ~, Jstat] = NSOLVE(@(U) GM.STATRESFUN(U, Fv*Prestress), ...
                                             U0, opts); 
        if eflag<=0
            error('No Prestress Convergence')
        end
    end
    Tstat = GM.NLTs.func(0, GM.NLTs.L*Ustat);
    save('statsolguess.mat', 'Ustat');
    
    %% Linearized Analysis
    [Vstat, Wstat] = eigs(Jstat, M, 10, 'SM');
    [Wstat, si] = sort(sqrt(diag(Wstat)));
    Vstat = Vstat(:, si);
    Vstat = Vstat./sqrt(diag(Vstat'*M*Vstat))';
    
    %% March
    Na = 10;
    As = logspace(min(Arange), max(Arange), Na);
    As = [-As(end:-1:1) As]';
    Eflags = zeros(1, 2*Na);

    UlC = zeros(GM.Ndofs+1, 2*Na);
    dUdalC = zeros(GM.Ndofs+1, 2*Na);

    ul0 = [Ustat+Vstat(:,mds(mdi))*As(Na+1); Wstat(mds(mdi))^2];
    opts = struct('reletol', 1e-6, 'Display', false, 'lsrch', 0);
    for ia=1:Na
        [UlC(:, Na+ia), ~, Eflags(Na+ia), ~, dRdUl] = ...
            NSOLVE(@(ul) GM.RQMRESFUN([ul; log10(As(Na+ia))], 1, ...
                                      Fv*Prestress, Ustat), ul0, opts);
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
        [UlC(:, Na+1-ia), ~, Eflags(Na+1-ia), ~, dRdUl] = ...
            NSOLVE(@(ul) GM.RQMRESFUN([ul; log10(-As(Na+1-ia))], 1, ...
                                      Fv*Prestress, Ustat), ul0, opts);
        
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5. Post Processing (Hermite Interpolation + SHBM) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ln = reshape([UlC(end,:); dUdalC(end,:)], 2*size(UlC,2), 1);

    Nq = 100;
    Qs = logspace(min(Qrange), max(Qrange), Nq)';

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
    Dfluxes = zeros(MESH.Nq^2*MESH.Ne*3, Nq);
    parfor (qi=1:Nq,8)
%     for qi=1:Nq        
        [~, Fnl] = GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)), tol, Nits);
        Dfluxes(:, qi) = mean((squeeze(Udot(:, qi, :))*GM.NLTs.L').*Fnl);
        
        Zts(qi) = GETFOURIERCOEFF(0, sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M +...
            squeeze(Udot(:, qi, :))*GM.C +...
            squeeze(Ut(:, qi, :))*GM.K + ...
            Fnl*GM.NLTs.Lf'),2))/(Qs(qi)^2*Lams(qi)^1.5);
        fprintf('%d\n', qi)
    end

    %% Create output structure & save to file
    
    save(sprintf('./%s_%d_m%d.mat', pref, nxi, mdi), 'Qs', ...
        'Phi', 'Lams', 'Zts', 'Wstat', 'Ustat', 'Xis', 'Wis', ...
        'Tstat', 'Dfluxes');
    outputs = load(sprintf('./%s_%d_m%d.mat', pref, nxi, mdi), 'Qs', ...
        'Phi', 'Lams', 'Zts', 'Wstat', 'Ustat', 'Xis', 'Wis', ...
        'Tstat', 'Dfluxes');   
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('=============================================\n')
    fprintf('Done %d\n', nxi);
    fprintf('=============================================\n')
    toc
end
