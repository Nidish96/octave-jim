clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/FEM/')

%% Setup System
M = [1 0;0 1];
K = [2 -1;-1 2];
L = [1 0];
% Lf = [0.5; 0];
% fnl = @(t, u, ud) deal(u.^3, 3*u.^2, zeros(size(u)));
fnl = @(t, u, ud) deal(0.5*u.^3, 1.5*u.^2, zeros(size(u)));
%% Setup Model
GM = MDOFGEN(M, K, M*0, eye(2));
% GM = GM.SETNLFUN(1+5, L, fnl, Lf);
GM = GM.SETNLFUN(1+3, L, fnl);

%% Linear Modes
[V, Wsp] = eig(K, M);
[Wsp,si] = sort(sqrt(diag(Wsp)));
V = V(:,si);
V = V./sqrt(diag(V'*M*V))';

%% Continuation
Amax = 2.5;
da = 0.1;

ul0 = [V(:,1)*10^Amax; Wsp(1)^2];
Dscale = [ones(GM.Ndofs, 1); Wsp(1)^2; 1.0];

Copt = struct('Nmax', 1000, 'Dscale', Dscale, 'dsmax', 0.45);
% [UlC, dUlC, Ss] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0), ul0, -10^Amax, 10^Amax, da, Copt);

ul0 = [V(:,1)*1e-2; Wsp(1)^2];
[UlC1, dUlC1, Ss] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0), ul0, 1e-2, 10^Amax, da, Copt);
ul0 = [-V(:,1)*1e-2; Wsp(1)^2];
[UlC2, dUlC2, Ss] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0), ul0, -1e-2, -10^Amax, da, Copt);

UlC = [UlC2(:, end:-1:1), UlC1(:, 1:end)];
dUlC = [dUlC2(:, end:-1:1), dUlC1(:, 1:end)];

%% Post Processing (Hermite interpolation + 1HB)
Ln = reshape([UlC(end-1,:); dUlC(end-1,:)], 2*size(UlC,2), 1);
As = UlC(end,:)';

Nq = 100;
% Qs = (10^Amax)*linspace(0.01, 1, Nq)';
Qs = logspace(-2, Amax, Nq)';

Nt = 2^7;
t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
qt = cos(t).*Qs';
tic
[Lt, Nint, dNint] = HERMINTERP(As, Ln, qt(:));
toc
Lt = reshape(Lt, Nt, Nq);

% Natural Frequency
Lams = sqrt(sum((GETFOURIERCOEFF(1, Lt.*qt)./Qs').^2))';
% sum(GETFOURIERCOEFF(1, Lt.*qt).*GETFOURIERCOEFF(1, Lt.*qt), 2)
qdot = -sin(t).*(sqrt(Lams).*Qs)';
qddot = -cos(t).*(Lams.*Qs)';

% Mode Shape
Un = reshape([permute(UlC(1:end-2,:), [3, 2, 1]); 
    permute(dUlC(1:end-2,:), [3, 2, 1])], [], GM.Ndofs);  % (2Npt,Ndofs)
Ut = reshape(Nint*Un, Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs)
Udot = reshape((dNint*Un).*qdot(:), Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs)
Uddot = reshape((dNint*Un).*qddot(:), Nt, Nq, GM.Ndofs);  % (Nt, Nq, Ndofs)

Uh = reshape(GETFOURIERCOEFF(1, reshape(Ut, Nt, Nq*GM.Ndofs)), 2, Nq, GM.Ndofs);
Phi = (squeeze(Uh(1,:,:)-1j*Uh(2,:,:))./Qs)';

% Damping
tic
Zts = zeros(Nq, 1);
fnl = @(t, u, ud) 0.5*u.^3;
GM.NLTs(1).func = fnl;
parfor (qi=1:Nq, 8)
%     size(sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M + squeeze(Udot(:, qi, :))*GM.C + squeeze(Ut(:, qi, :))*GM.K + GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)))),2))
    Zts(qi) = GETFOURIERCOEFF(0, sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M +...
        squeeze(Udot(:, qi, :))*GM.C +...
        squeeze(Ut(:, qi, :))*GM.K + GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)))),2))/(Qs(qi)^2*Lams(qi)^1.5);
end
toc

%% Save
save('./DATA/Stifnl_RQNM.mat', 'Qs', 'Zts', 'Lams', 'Phi', 'UlC', 'dUlC', 'GM', 'Nint', 'dNint');

%% Load
load('./DATA/Stifnl_RQNM.mat', 'Qs', 'Zts', 'Lams', 'Phi')
load('./DATA/Stifnl_EPMC.mat', 'UwxC');
%% Plot
figure(1)
clf()
hold on
semilogx(Qs, sqrt(Lams), '.-'); hold on
plot(10.^UwxC(end,:), UwxC(end-2,:), '.-'); hold on
set(gca, 'XScale', 'log')

xlabel('Modal Amplitude')
ylabel('Frequency (rad/s)')

figure(2)
clf()
hold on
plot(Qs, Zts, '.-')

% set(gca, 'YScale', 'log')
xlabel('Modal Amplitude')
ylabel('Damping Factor')

figure(3)
clf()
br=bar3(1-abs((Phi'*GM.M*Phi).^2./(diag(Phi'*GM.M*Phi).*diag(Phi'*GM.M*Phi)')));

for k=1:size(br,2); br(k).CData = br(k).ZData; br(k).FaceColor = 'interp'; br(k).EdgeColor = 'none'; end