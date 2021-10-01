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

[V, Wsp] = eig(Kb, Mb);
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
cfg = 1;

if cfg==1
    % Shows 1:3 mode coupling (cfg=1)
    kt  = 6e6;
    muN = 7.5e5;
elseif cfg==2
    % Doesn't show mode coupling (cfg=2)
    kt  = 2.5e6;
    muN = 7.5e5;
else
    error('unknown cfg')
end


fnl = @(t, u, varargin) JENKFORCE(t, u, kt, muN, varargin{:});
GM = GM.SETNLFUN(2+3, Lb(end,:), fnl);

%% Continuation

Amax = 2;
da = 0.01;

mds = [1 2];
Copt = struct('Nmax', 10000, 'dsmax', 0.2);
UlC = cell(size(mds));
dUlC = cell(size(mds));
Ss = cell(size(mds));
for mi=mds
    ul0 = [V(:,mi)*10^Amax; Wsp(mi)^2];
    Dscale = [ones(GM.Ndofs, 1); Wsp(mi)^2; 1.0];
    Copt.Dscale = Dscale;
    
    [UlC{mi}, dUlC{mi}, Ss{mi}] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0), ul0, -10^Amax, 10^Amax, da, Copt);
end
%% Post-Processing for Dynamic Modal Properties
Nq = 100;
% Qs = (10^Amax)*linspace(0.01, 1, Nq)';
Qs = logspace(-2, Amax, Nq)';

Nt = 2^7;
t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
qt = cos(t).*Qs';

Lams = zeros(Nq, length(mds));
Cs = zeros(Nq, length(mds));
Omegas = zeros(Nq, length(mds));
for mi=mds
    Ln{mi} = reshape([UlC{mi}(end-1,:); dUlC{mi}(end-1,:)], 2*size(UlC{mi},2), 1);
    As{mi} = UlC{mi}(end,:)';

    tic
    [Lt{mi}, Nint{mi}, dNint{mi}] = HERMINTERP(As{mi}, Ln{mi}, qt(:));
    toc
    Nint{mi} = sparse(Nint{mi});  % Interpolant
    dNint{mi} = sparse(dNint{mi});  % Gradient Interpolant
    
    Lt{mi} = reshape(Lt{mi}, Nt, Nq);
    
    % Natural Frequency
	FH = GETFOURIERCOEFF(1, Lt{mi}.*qt);  % Forcing harmonics
    Omegas(:, mi) = sqrt(FH(1,:)./Qs');  % Cosine force harmonic/Q: LS fit in freq. domain
    qdot = -sin(t).*(Omegas(:, mi).*Qs)';
    qddot = -cos(t).*(Omegas(:, mi).^2.*Qs)';
    
    % Mode Shape
    Un{mi} = reshape([permute(UlC{mi}(1:end-2,:), [3, 2, 1]); 
        permute(dUlC{mi}(1:end-2,:), [3, 2, 1])], [], GM.Ndofs);  % (2Npt,Ndofs)
    Ut{mi} = reshape(Nint{mi}*Un{mi}, Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs): Time domain deflection
    Udot{mi} = reshape((dNint{mi}*Un{mi}).*qdot(:), Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs): Time domain velocity
    Uddot{mi} = reshape((dNint{mi}*Un{mi}).*qddot(:), Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs): Time domain velocity
    dUdqt{mi} = reshape(dNint{mi}*Un{mi}, Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs): Time domain deflection
    
%     Phi{mi} = reshape(GETFOURIERCOEFF(1, reshape(Ut{mi}, Nt, Nq*GM.Ndofs)), 2, Nq, GM.Ndofs)./Qs';
%     Phi{mi} = reshape(GETFOURIERCOEFF(0, reshape(dUdqt{mi}, Nt, Nq*GM.Ndofs)), 1, Nq, GM.Ndofs);
    Phi{mi} = reshape(GETFOURIERCOEFF(0, reshape(dUdqt{mi}, Nt, Nq*GM.Ndofs)), Nq, GM.Ndofs);  % (Nq, Ndofs): Real Mode shape
    
    % Single mode dissipation factor
    tol = 1e-6;
    for qi=1:Nq
        force_nl = GM.NLEVAL(t, squeeze(Ut{mi}(:, qi, :)), squeeze(Udot{mi}(:, qi, :)), tol) + ...
            squeeze(Udot{mi}(:, qi, :))*GM.C' + squeeze(Ut{mi}(:, qi, :))*GM.K' + squeeze(Uddot{mi}(:, qi, :))*GM.M';
%         Cs(qi, mi) = 2*GETFOURIERCOEFF(0, sum(squeeze(Udot{mi}(:, qi, :)).*force_nl, 2))/(Qs(qi)^2*Omegas(qi,mi)^2);  % From energy dissipated
        Cs(qi, mi) = [0;-Omegas(qi,mi)*Qs(qi)]\GETFOURIERCOEFF(1,force_nl*Phi{mi}(qi,:)');  % Regression
    end
end

% Dissipation
Ncyc = 5;
fsamp = 2^(ceil(log2(max(Omegas(:))/2/pi))+4);
T = (0:1/fsamp:max(2*pi./Omegas(:))*Ncyc)';
force_t = zeros(length(T), GM.Ndofs);
Ns = fix(max(2*pi./Omegas(:))*fsamp*2);  % Start for regression

if length(mds)==2
    Cnl = zeros(3, Nq, Nq);
    Knl = zeros(3, Nq, Nq);
    Q1s = zeros(1, Nq, Nq);
    Q2s = zeros(1, Nq, Nq);
    for q1i=1:Nq
        q1t = cos(Omegas(q1i, 1)*T)*Qs(q1i);
        q1dot = -Omegas(q1i, 1)*sin(Omegas(q1i, 1)*T)*Qs(q1i);
        q1ddot = -Omegas(q1i, 1)^2*cos(Omegas(q1i, 1)*T)*Qs(q1i);
        Q1s(:, q1i, :) = Qs(q1i);
        for q2i=1:Nq
            q2t = cos(Omegas(q2i, 2)*T)*Qs(q2i);
            q2dot = -Omegas(q2i, 2)*sin(Omegas(q2i, 2)*T)*Qs(q2i);
            q2ddot = -Omegas(q2i, 2)^2*cos(Omegas(q2i, 2)*T)*Qs(q2i);
            Q2s(:, :, q2i) = Qs(q2i);
            
            ut = Phi{1}(q1i,:).*q1t + Phi{2}(q2i,:).*q2t;
            udot = Phi{1}(q1i, :).*q1dot + Phi{2}(q2i, :).*q2dot;
            
            % Initialize slider states
            GM.NLTs.up = 0;
            GM.NLTs.fp = 0;
            for ti=1:length(T)
                [force_t(ti,:), ~, ~, GM] = GM.NLFORCE(T(ti), ut(ti,:)', udot(ti,:)', T(ti)-T(2));  % Nonlinear force
                force_t(ti, :) = force_t(ti, :) + ...
                    udot(ti, :)*GM.C' + ut(ti, :)*GM.K'; % Adding Linear forcing dissipation
            end
            
            % Regression to fit c11, c12, c22 (symmetric modal C)
            fmodal_t = ([Phi{1}(q1i, :); Phi{2}(q2i, :)]*force_t')';  % (Nt, 2)
            
            % only fit cnl
%             fmodal_t = fmodal_t - [q1t q2t]*diag([Omegas(q1i, 1)^2; Omegas(q2i, 2)^2]);  % Removing "conservative part"            
% %             Cnl(:, q1i, q2i) = [q1dot q2dot zeros(length(T),1);zeros(length(T),1) q1dot q2dot]\[fmodal_t(:,1); fmodal_t(:,2)]; % [c11; c12; c22]
%             Cnl(:, q1i, q2i) = [q1dot(Ns:end) q2dot(Ns:end) zeros(length(T)-Ns+1,1);...
%                 zeros(length(T)-Ns+1,1) q1dot(Ns:end) q2dot(Ns:end)]\[fmodal_t(Ns:end,1); fmodal_t(Ns:end,2)]; % [c11; c12; c22]
            
            % fit knl, cnl
            reg = [[q1t(Ns:end) q2t(Ns:end) zeros(length(T)-Ns+1,1);...
                zeros(length(T)-Ns+1,1) q1t(Ns:end) q2t(Ns:end)], ...
                [q1dot(Ns:end) q2dot(Ns:end) zeros(length(T)-Ns+1,1);...
                zeros(length(T)-Ns+1,1) q1dot(Ns:end) q2dot(Ns:end)]]\[fmodal_t(Ns:end,1); fmodal_t(Ns:end,2)];
            Knl(:, q1i, q2i) = reg(1:3);
            Cnl(:, q1i, q2i) = reg(4:6);
            
            fprintf('(%d, %d)\n', q1i, q2i);
        end
    end
else
    error('Code only written for 2 modes under consideration')
end

%% Assess Fits
q1i = 1;
q2i = 50;

KK = Knl(:, q1i, q2i);
KK = [KK(1) KK(2);KK(2) KK(3)];
CC = Cnl(:, q1i, q2i);
CC = [CC(1) CC(2);CC(2) CC(3)];

q1t    = cos(Omegas(q1i, 1)*T)*Qs(q1i);
q1dot  = -Omegas(q1i, 1)*sin(Omegas(q1i, 1)*T)*Qs(q1i);
q1ddot = -Omegas(q1i, 1)^2*cos(Omegas(q1i, 1)*T)*Qs(q1i);

q2t    = cos(Omegas(q2i, 2)*T)*Qs(q2i);
q2dot  = -Omegas(q2i, 2)*sin(Omegas(q2i, 2)*T)*Qs(q2i);
q2ddot = -Omegas(q2i, 2)^2*cos(Omegas(q2i, 2)*T)*Qs(q2i);

ut = Phi{1}(q1i,:).*q1t + Phi{2}(q2i,:).*q2t;
udot = Phi{1}(q1i, :).*q1dot + Phi{2}(q2i, :).*q2dot;
            
% Initialize slider states
GM.NLTs.up = 0;
GM.NLTs.fp = 0;
for ti=1:length(T)
    [force_t(ti,:), ~, ~, GM] = GM.NLFORCE(T(ti), ut(ti,:)', udot(ti,:)', T(ti)-T(2));  % Nonlinear force
    force_t(ti, :) = force_t(ti, :) + udot(ti, :)*GM.C' + ut(ti, :)*GM.K'; % Adding Linear forcing dissipation
end
            
% Regression to fit c11, c12, c22 (symmetric modal C)
fmodal_t = ([Phi{1}(q1i, :); Phi{2}(q2i, :)]*force_t')';  % (Nt, 2)

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',16)
figure(10)
clf()
predfm_t = [q1t q2t]*KK+[q1dot q2dot]*CC;
for si=1:2
    subplot(2,1, si)
%     plot(T, fmodal_t(:, si)); hold on
%     plot(T, fmodal_t(:,si)-predfm_t(:,si), 'k-')

    plot(T, fmodal_t(:, si), 'k.'); hold on
    plot(T, predfm_t(:, si), 'b-')
    
    ylabel(sprintf('$f^{(m,%d)}$', si))
end
xlabel('Time (s)')
% print(sprintf('./FIGS/fitdepict_%d,%d.eps', q1i,q2i), '-depsc')
%% Save
% save(sprintf('./DATA/Jenkins_EPMC_cfg%d.mat',cfg), 'UwxC', 'dUwxC', 'h', 'Nhc', 'mds');
save(sprintf('./DATA/Jenkins_RQNM_cfg%d.mat',cfg), 'UlC', 'dUlC', 'As', 'Phi', 'Cs', 'Omegas', 'Q1s', 'Q2s', 'Knl', 'Cnl', 'GM', 'Qs')

%% plot
cfg = 1;
load(sprintf('./DATA/Jenkins_RQNM_cfg%d.mat',cfg), 'UlC', 'dUlC', 'As', 'Phi', 'Cs', 'Omegas', 'Q1s', 'Q2s', 'Knl', 'Cnl', 'GM', 'Qs')
% [Q2s, Q1s] = meshgrid(Qs, Qs);

figure(cfg*100)
clf();

subplot(4,2, 1)
surf(squeeze(Q1s), squeeze(Q2s), squeeze(Knl(1, :, :)), 'EdgeColor', 'None');
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
% xlabel('Q1')
ylabel('Q2')
colorbar
title('K(1,1)')
set(gca, 'view', [0 90])
caxis([6.5 7.3]*1e7)

subplot(4,2, 3)
surf(squeeze(Q1s), squeeze(Q2s), squeeze(Knl(2, :, :)), 'EdgeColor', 'None');
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
% xlabel('Q1')
ylabel('Q2')
colorbar
title('K(1,2)')
set(gca, 'view', [0 90])
caxis([-5.7 1.5]*1e6)

subplot(4,2, 5)
surf(squeeze(Q1s), squeeze(Q2s), squeeze(Knl(3, :, :)), 'EdgeColor', 'None');
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
xlabel('Q1')
ylabel('Q2')
colorbar
title('K(2,2)')
set(gca, 'view', [0 90])
caxis([6 6.5]*1e8)

subplot(4,2, 2)
surf(squeeze(Q1s), squeeze(Q2s), squeeze(Cnl(1, :, :)), 'EdgeColor', 'None');
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
% xlabel('Q1')
% ylabel('Q2')
colorbar
title('C(1,1)')
set(gca, 'view', [0 90])
caxis([65 425])

subplot(4,2, 4)
surf(squeeze(Q1s), squeeze(Q2s), squeeze(Cnl(2, :, :)), 'EdgeColor', 'None');
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
% xlabel('Q1')
% ylabel('Q2')
colorbar
title('C(1,2)')
set(gca, 'view', [0 90])
caxis([-1.5 212])

subplot(4,2, 6)
surf(squeeze(Q1s), squeeze(Q2s), squeeze(Cnl(3, :, :)), 'EdgeColor', 'None');
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
xlabel('Q1')
% ylabel('Q2')
colorbar
title('C(2,2)')
set(gca, 'view', [0 90])
caxis([-600 600])

subplot(4,4, 14:15)
surf(squeeze(Q1s), squeeze(Q2s), Phi{1}*GM.M*Phi{2}', 'EdgeColor', 'None')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
xlabel('Q1')
ylabel('Q2')
colorbar
title('m(1,2)')
set(gca, 'view', [0 90])
caxis(1e-2*[-1 1])
return

%%
figure(300)
surf(squeeze(Q1s), squeeze(Q2s), squeeze(Cnl(1, :, :).*(Q1s.*Omegas(:,1)).^2/2) + squeeze(Cnl(2, :, :)).*(Q2s.*Omegas(:,2)').^2/2)

set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
xlabel('Q1')
ylabel('Q2')