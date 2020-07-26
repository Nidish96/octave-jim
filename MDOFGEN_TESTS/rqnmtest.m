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

kc = 1e6;
fnl = @(t,u,ud) deal(kc*u.^3, 3*kc*u.^2, zeros(size(u)));
GM = GM.SETNLFUN(1+3, Lb(end,:), fnl);

%% Continuation

Amax = 0.5;
da = 0.1;

ul0 = [V(:,1)*10^Amax; Wsp(1)^2];
Dscale = [ones(GM.Ndofs, 1); Wsp(1)^2; 1.0];

Copt = struct('Nmax', 100, 'Dscale', Dscale);
[UlC, dUlC, Ss] = CONTINUE(@(ulq) GM.RQMRESFUN(ulq,0), ul0, -10^Amax, 10^Amax, da, Copt);

%% Post Processing (Hermite interpolation + 1HB)
Ln = reshape([UlC(end-1,:); dUlC(end-1,:)], 2*size(UlC,2), 1);
As = UlC(end,:)';

Nq = 100;
Qs = floor(10^Amax)*linspace(0.01, 1, Nq)';

Nt = 2^7;
t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
qt = cos(t).*Qs';
[Lt, Nint, dNint] = HERMINTERP(As, Ln, qt(:));
Lt = reshape(Lt, Nt, Nq);

% Natural Frequency
Lams = sqrt(sum((GETFOURIERCOEFF(1, Lt.*qt)./Qs').^2))';
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

Zts = zeros(Nq, 1);
for qi=1:Nq
%     size(sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M + squeeze(Udot(:, qi, :))*GM.C + squeeze(Ut(:, qi, :))*GM.K + GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)))),2))
    Zts(qi) = GETFOURIERCOEFF(0, sum(squeeze(Udot(:, qi, :)).*(squeeze(Uddot(:, qi, :))*GM.M +...
        squeeze(Udot(:, qi, :))*GM.C +...
        squeeze(Ut(:, qi, :))*GM.K + GM.NLEVAL(t, squeeze(Ut(:, qi,:)), squeeze(Udot(:, qi,:)))),2))/(Qs(qi)^2*Lams(qi)^1.5);
end
disp('done')

%%
figure(1)
% clf()
% plot(sqrt(UlC(end-1,:)), UlC(end,:), 'o-'); hold on
plot(sqrt(Lams), Qs, '.-')

xlabel('Frequency(Hz)')
ylabel('Amplitude')

figure(2)
plot(Qs, Zts, '.-')

xlabel('Amplitude')
ylabel('Damping Factor')


figure(3)
clf()
br=bar3(1-abs((Phi'*GM.M*Phi).^2./(diag(Phi'*GM.M*Phi).*diag(Phi'*GM.M*Phi)')));

for k=1:size(br,2); br(k).CData = br(k).ZData; br(k).FaceColor = 'interp'; br(k).EdgeColor = 'none'; end