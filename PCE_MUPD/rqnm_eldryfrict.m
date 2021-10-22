clc
clear all

addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/CONTACTMODELS/')
addpath('../ROUTINES/QUASISTATIC/')
addpath('../ROUTINES/SOLVERS/')

model = 'BRB';

load(sprintf('../MODELS/%s/MATRICES_NR.mat', model), 'K', 'M', 'R', 'L', 'Fv', 'TFM');

%% Load Mesh
Nds = dlmread(sprintf('../MODELS/%s/Nodes.dat', model));
Quad = dlmread(sprintf('../MODELS/%s/Elements.dat', model));

Nq = 2;
MESH = MESH2D(Nds, 3, [], Quad, Nq);

%% Create Object
GM = MDOFGEN(M, K, zeros(size(M)), L);

Kt = [1e12; 1e12; 0];
kn = 1e12;
mu = 0.25;
gap = 0;

fnl = @(t, u, varargin) ELDRYFRICT(t, u, Kt, kn, mu, gap, varargin{:});
GM = GM.SETNLFUN(2+5, ...
    kron(MESH.Qm, eye(3))*L(1:(MESH.Nn*3), :), ...
    fnl, ...
    L(1:(MESH.Nn*3), :)'*kron(MESH.Tm, eye(3)));

%% Fully Engaged Stiffness
K0 = L(1:MESH.Nn*3,:)'*kron(MESH.Tm, eye(3))*kron(eye(MESH.Ne*MESH.Nq^2), [Kt(1) Kt(3) 0; Kt(3) Kt(2) 0; 0 0 kn])*kron(MESH.Qm, eye(3))*L(1:MESH.Nn*3,:);
% K0 = (K0+K0')/2;

% K1 = zeros(size(L,1));
% K1(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Kt(1);
% K1(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Kt(2);
% K1(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*kn;
% K1 = L'*K1*L;

%% Static Prestress
opts = struct('reletol', 1e-6, 'Display', true, 'lsrch', 0);

Prestress = 12e3;
U0 = (K+K0)\(Fv*Prestress);

[Ustat, ~, ~, ~, Jstat] = NSOLVE(@(U) GM.STATRESFUN(U, Fv*Prestress), U0, opts);
% U0 = Ustat;

% %% 
% % Tstat = GM.NLTs.func(0, GM.NLTs.L*Ustat);
% Tstat = GM.NLTs.func(0, L(1:MESH.Nn*3,:)*Ustat);
% 
% figure(1)
% clf()
% for i=1:3
%     subplot(3,1,i)
%     MESH.SHOWFIELD2D(Tstat(i:3:end))
%     xx=colorbar('southoutside');
%     xlabel(xx, sprintf('Traction %d', i))
%     axis equal
% end
%     

%% Linearized Analysis
[Vstat, Wstat] = eigs(Jstat, M, 10, 'SM');
[Wstat, si] = sort(sqrt(diag(Wstat)));
Vstat = Vstat(:, si);
Vstat = Vstat./sqrt(diag(Vstat'*M*Vstat))';

% %% RQNM Continuation
% Amax = -3;
% da = 0.10;
% 
% Astart = -6;
% 
% ul0 = [Vstat(:,1)*10^Amax; Wstat(1)^2];
% % [R, dR] = GM.RQMRESFUN([ul0; 10^Amax],0);
% 
% ul0 = [Ustat+Vstat(:,1)*10^Astart; Wstat(1)^2];
% opts = struct('reletol', 1e-6, 'Display', true, 'lsrch', 0);
% % Usol = NSOLVE(@(ul) GM.RQMRESFUN([ul; 10^Astart],0,Fv*Prestress,Ustat), ul0, opts);
% % Usol = NSOLVE(@(ul) GM.RQMRESFUN([ul; Astart],1,Fv*Prestress,Ustat), ul0, opts);
% 
% disp([Wstat(1) sqrt(Usol(end))]/2/pi)
% % opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% % Usol = fsolve(@(ul) GM.RQMRESFUN([ul; 10^da],0,Prestress*Fv,Ustat), ul0, opts);
% % Dscale = [ones(GM.Ndofs, 1); Wstat(1)^2; 1.0];
% Dscale = [ones(size(Usol));1];
% % RQMRESFUN
% 
% ul0 = Usol;
% Copt = struct('Nmax', 10000, 'Dscale', Dscale, 'dsmax', 0.2);
% % [UlC, dUlC, Ss] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0,Fv*Prestress,Ustat), ul0, 10^Astart, 10^Amax, da, Copt);
% da = 0.5;
% % UlC = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,1,Fv*Prestress,Ustat), ul0, Astart, Amax, da, Copt);

%% March
Na = 5;
As = logspace(-6, -3, Na);
As = [-As(end:-1:1) As];
UlC = zeros(GM.Ndofs+2, 2*Na);

tic
ul0 = [Ustat+Vstat(:,1)*As(Na+1); Wstat(1)^2];
opts = struct('reletol', 1e-6, 'Display', true, 'lsrch', 0);
% UlC(:,1) = NSOLVE(@(ul) GM.RQMRESFUN([ul; Astart],1,Fv*Prestress,Ustat), ul0, opts);
for ia=1:Na
%     UlC(:, Na+ia) = [NSOLVE(@(ul) GM.RQMRESFUN([ul; log10(As(Na+ia))], 1, Fv*Prestress, Ustat), ul0, opts); As(Na+ia)];
    UlC(:, Na+ia) = [NSOLVE(@(ul) GM.RQMRESFUN([ul; As(Na+ia)], 0, Fv*Prestress, Ustat), ul0, opts); As(Na+ia)];
    if ia<Na
        ul0 = UlC(1:end-1, Na+ia);
        ul0(1:end-1) = Ustat+(ul0(1:end-1)-Ustat)*As(Na+ia+1)/As(Na+ia);
    end
    
    fprintf('%d, ', ia);
end
fprintf('\n');

ul0 = [Ustat+Vstat(:,1)*As(Na+1-1); Wstat(1)^2];
for ia=1:Na
    UlC(:, Na+1-ia) = [NSOLVE(@(ul) GM.RQMRESFUN([ul; -(As(Na+1-ia))], 0, Fv*Prestress, Ustat), ul0, opts); As(Na+1-ia)];
    if ia<Na
        ul0 = UlC(1:end-1, Na+1-ia);
        ul0(1:end-1) = Ustat+(ul0(1:end-1)-Ustat)*As(Na+1-ia-1)/As(Na+1-ia);
    end
    
    fprintf('%d, ', ia);
end
toc

%% Plot Backbone
