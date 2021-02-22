clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/FEM/')

%% Setup system (1D Bar)
ell = 1.0;
A = pi*1e-4;
E = 210e9;
rho = 7800;
Ne = 8;

Le = ell/Ne;
Me = 1/3*[2 1;1 2]*Le/2*rho*A;
Ke = 1/2*[1 -1;-1 1]*2/Le*E*A;

M = zeros(Ne+1);
K = zeros(Ne+1);
for e=1:Ne
    M(e:e+1,e:e+1) = M(e:e+1,e:e+1) + Me;
    K(e:e+1,e:e+1) = K(e:e+1,e:e+1) + Ke;
end
Lb = eye(Ne+1);
Lb(:, 1) = [];  % Boundary Condition (DOF 1 fixed)

Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;

%% Damping Factors
Zetas = [0.4; 0.4]*1e-2;

[V, Wsp] = eig(Kb, Mb);  % Slipped modes 
[Wsp, si] = sort(sqrt(diag(Wsp)));
V = V(:, si);
V = V./diag(sqrt(V'*Mb*V))';

if Ne>2
    ab = [1./(2*Wsp(1:length(Zetas))) Wsp(1:length(Zetas))/2]\Zetas;
    Cb = ab(1)*Mb + ab(2)*Kb;
else
    Cb = 2*Zetas(1)/Wsp(1)*Kb;
end

%% Setup model
GM = MDOFGEN(Mb, Kb, Cb, Lb);

% kc = 1e6;
% fnl = @(t,u,ud) deal(kc*u.^3, 3*kc*u.^2, zeros(size(u)));
kt  = 7.5e6;
muN = 7.5e5;
fnl = @(t, u, varargin) JENKFORCE(t, u, kt, muN, varargin{:});
GM = GM.SETNLFUN(2+3, Lb(end,:), fnl);

%% Continuation for RQNM backbone
Amax = 2;
da = 5;

% % I'm starting at -10^Amax and going to 10^Amax here (might not work always)
% ul0 = [V(:,1)*10^Amax; Wsp(1)^2];
% Dscale = [abs(ul0); 1.0];
% Copt = struct('Nmax', 10000, 'Dscale', Dscale, 'dsmax', 20, 'DynDscale', 1);
% Copt = struct('Nmax', 10000, 'dsmax', 20);
% [UlC, dUlC, Ss] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0), ul0, -10^Amax, 10^Amax, da, Copt);

% % For cases when the first one doesn't converge/struggles around the origin, do this:
startamp = 1e-2;

ul0 = [V(:, 1)*startamp; Wsp(1)^2];
Dscale = [abs(ul0); 1.0];
Copt = struct('Nmax', 10000, 'Dscale', Dscale, 'dsmax', 20, 'DynDscale', 1);
[UlC1, dUlC1, Ss] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0), ul0, startamp, 10^Amax, da, Copt);
ul0 = [-V(:, 1)*startamp; Wsp(1)^2];
[UlC2, dUlC2, Ss] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0), ul0, -startamp, -10^Amax, da, Copt);

UlC = [UlC2(:, end:-1:1), UlC1(:, 1:end)];
dUlC = [dUlC2(:, end:-1:1), dUlC1(:, 1:end)];

%% Post Processing (Hermite interpolation + 1HB)
Ln = reshape([UlC(end-1,:); dUlC(end-1,:)], 2*size(UlC,2), 1);
As = UlC(end,:)';
% Mode Shape
Un = reshape([permute(UlC(1:end-2,:), [3, 2, 1]); 
    permute(dUlC(1:end-2,:), [3, 2, 1])], [], GM.Ndofs);  % (2Npt,Ndofs)

Nq = 100;
% Qs = (10^Amax)*linspace(0.01, 1, Nq)';
Qs = logspace(-2, Amax, Nq)';

Nt = 2^7;  % Time points for AFT
t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
qw20 = [0; Wsp(1)^2];

opt = optimset('Jacobian', 'on', 'Display', 'off');

% Nonlinear harmonic force in subspace is fm(t)=\lambda(q(t))*q(t)
% Flow rule is qdd+fm(t)=0

% "Conservative Properties"
qw2s = zeros(3, Nq);

dUt = zeros(Nt, Nq, GM.Ndofs);
Ut = zeros(Nt, Nq, GM.Ndofs);
Udot = zeros(Nt, Nq, GM.Ndofs);
Uddot = zeros(Nt, Nq, GM.Ndofs);
for qi=1:Nq
    qw2s(1:end-1, qi) = fsolve(@(qw2) RQHBRESFUN([qw2; Qs(qi)], @(t, q) HERMINTERP(As, Ln, q(:)).*q(:), ...
        @(t, q) HERMINTERP(As, Ln, q(:))+interp1(As, Ln(2:2:end), q).*q(:), Nt), qw20, opt);
    qw2s(end, qi) = Qs(qi);
    qw20 = qw2s(1:end-1, qi);
    
    % Mode Shape
    qt = qw2s(1, qi)+qw2s(3, qi)*cos(t);
    [~, Nint, dNint] = HERMINTERP(As, Ln, qt);
    
    Ut(:, qi, :) = Nint*Un;
    dUt(:, qi, :) = dNint*Un;
    Udot(:, qi, :) = (dNint*Un).*(-sqrt(qw2s(2,qi))*qw2s(3,qi)*sin(t));
%     Uddot(:, i, :) = (dNint*Un).*qddot(:, i);
    
    fprintf('%d/%d\n', qi, Nq)
end
Lams = qw2s(end-1, :);
qdot = -sin(t).*(sqrt(Lams(:)).*Qs(:))';
qddot = -cos(t).*(Lams(:).*Qs(:))';

% U/Q formulation
Uh = reshape(GETFOURIERCOEFF(1, reshape(Ut, Nt, Nq*GM.Ndofs)), 2, Nq, GM.Ndofs);
Phi1 = (squeeze(Uh(1,:,:)-1j*Uh(2,:,:))./Qs)';

% dU/dq formulation
dUh = squeeze(reshape(GETFOURIERCOEFF(0, reshape(dUt, Nt, Nq*GM.Ndofs)), 1, Nq, GM.Ndofs))';
Phi2 = dUh;

% Normalize (not really necessary)
Phi1 = Phi1./sqrt(diag(Phi1'*GM.M*Phi1)');
Phi2 = Phi2./sqrt(diag(Phi2'*GM.M*Phi2)');

% Damping
tol = 1e-6;
tic
Zt_mf1 = zeros(Nq, 1);
Zt_mf2 = zeros(Nq, 1);
Zt_ds = zeros(Nq, 1);
% parfor (qi=1:Nq,8)  % Can be task parallelized
for qi=1:Nq
    Fint = (squeeze(Udot(:, qi, :))*GM.C +...
            squeeze(Ut(:, qi, :))*GM.K + ...
            GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)), tol));
    
% Modal Force Concept    
% \int \dot{u}^T f dt = 
% \int \left({\frac{\partial u}{\partial q}}^T f\right) \dot{q} dt
% \approxeq \int \left(\phi^T f\right) \dot{q}
% \implies f_{modal} = {\frac{\partial u}{\partial q}}^T f \approxeq \phi^T f

% 1. U/Q concept
    fmodal = Fint*real(Phi1(:,qi));
    C = GETFOURIERCOEFF(1, fmodal); 
    C = C(2);  % Choose sine Harmonic
    
    Zt_mf1(qi) = C/(-2*Qs(qi)*Lams(qi));
% 2. dU/dQ concept
    fmodal = Fint*real(Phi2(:,qi));
    C = GETFOURIERCOEFF(1, fmodal); 
    C = C(2);  % Choose sine Harmonic
    
    Zt_mf2(qi) = C/(-2*Qs(qi)*Lams(qi));    

% Total Dissipation Concept
    Diss = mean(sum(squeeze(Udot(:, qi, :)).*Fint, 2));  % zero harmonic of f(t)*udot(t)
    
    Zt_ds(qi) = Diss/(Qs(qi)^2*Lams(qi)^1.5);
    
    fprintf('%d/%d\n', qi, Nq)
end
toc

%% Plotting

figure(1)
clf()
semilogx(Qs, sqrt(Lams), '.-')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (rad/s)')

figure(2)
clf()
semilogx(Qs, Zt_ds, 'k.-'); hold on
semilogx(Qs, Zt_mf1, 'o-')
semilogx(Qs, Zt_mf2, '+-')

legend('Total Dissipation Concept', 'Modal Force Concept: U/Q Mode shape', ...
    'Modal Force Concept: dU/dQ Mode shape', 'Location', 'southeast')
xlabel('Modal Amplitude')
ylabel('Damping Factor')

figure(3)
clf()
semilogx(Qs, Zt_ds, 'k.-'); hold on
semilogx(Qs, Zt_mf1, 'o-')
semilogx(Qs, Zt_mf2, '+-')

xlim([0.16 0.35])
legend('Total Dissipation Concept', 'Modal Force Concept: U/Q Mode shape', ...
    'Modal Force Concept: dU/dQ Mode shape', 'Location', 'southeast')
xlabel('Modal Amplitude')
ylabel('Damping Factor')