clc
clear all
addpath('../../ROUTINES/')
addpath('../../ROUTINES/SOLVERS/')
addpath('../../ROUTINES/HARMONIC/')

%% Parameters
Ne = 10;
rho = 7800;
E = 2e11;
A = 0.01;
L = 1.0;

Me = rho*A*L/(6*Ne)*[2 1;1 2];
Ke = A*E*Ne/L*[1 -1;-1 1];
M = zeros(Ne+1);
K = zeros(Ne+1);
for e=1:Ne
    M(e:e+1, e:e+1) = M(e:e+1, e:e+1) + Me;
    K(e:e+1, e:e+1) = K(e:e+1, e:e+1) + Ke;
end
Lb = eye(Ne+1);
Lb(:, 1) = [];

Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;

%% Set Damping (Rayleigh Prop)
[Vb, Ws] = eig(Kb, Mb);
[Ws, si] = sort(sqrt(diag(Ws)));

Vb = Vb(:, si);

zts = [0.1; 0.2]*1e-2;
ab = [1./(2*Ws(1:2)) Ws(1:2)/2]\zts;

Cb = ab(1)*Mb + ab(2)*Kb;

%% model
MDL = MDOFGEN(Mb, Kb, Cb, L);

%% Excitation
bw = 1000;

% IMPULSE
bw = 1000;
famp = 100;
type = 'IMP';
fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw));
FEX = @(t) Lb(end,:)'*fex(t);

%% March
fsamp = 2^16;
T0 = 0;  T1 = 1200/fsamp;  dt = 1./fsamp;

opts = struct('Display', 'progress');

[Trk, Urk, Udrk, Uddrk] = MDL.RKGENMARCH(T0, T1, dt, zeros(MDL.Ndofs,1), ...
    zeros(MDL.Ndofs,1), FEX, opts);
[Thh, Uhh, Udhh, Uddhh] = MDL.HHTAMARCH(T0, T1, dt, zeros(MDL.Ndofs,1), ...
    zeros(MDL.Ndofs,1), FEX, opts);

%% 
func = @(t, y) [zeros(Ne), eye(Ne); -Mb\Kb, -Mb\Cb]*y + [zeros(Ne,1); Mb\Lb(end,:)']*fex(t);

[To, yo] = ode45(func, T0:dt:T1, zeros(MDL.Ndofs*2,1));

%% Analytical
[~, Tha, Xha] = lsim(ss([zeros(Ne), eye(Ne); -Mb\Kb, -Mb\Cb], [zeros(Ne,1); Mb\Lb(end,:)'], [Lb(end,:) Lb(end,:)*0], 0), ...
    fex(Thh), Thh);
% [~, Thr, Xhr] = lsim(ss([zeros(Ne), eye(Ne); -Mb\Kb, -Mb\Cb], [zeros(Ne,1); Mb\Lb(end,:)'], [Lb(end,:) Lb(end,:)*0], 0), ...
%     fex(Trk), Trk);
%% Plots
figure(1)
clf()

plot(Trk, Lb(end,:)*Urk, '.-'); hold on
plot(Thh, Lb(end,:)*Uhh, '+-'); hold on
plot(Tha, Lb(end,:)*Xha(:, 1:Ne)', 'o-'); hold on
% plot(Thr, Lb(end,:)*Xhr(:, 1:Ne)', 'o-'); hold on
% plot(To, Lb(end,:)*yo(:, 1:Ne)', 'o')

xlabel('Time (s)')
ylabel('Displacement (m)')

yyaxis right
plot(Trk, fex(Trk), 'k--')

ylabel('Excitation')

figure(2)
clf()
Xh = Lb(end,:)*Xha(:, 1:Ne)';
Xr = interp1(Tha, Lb(end,:)*Xha(:, 1:Ne)', Trk);

plot(Thh, Lb(end,:)*Uhh-Xh, '.-'); hold on
plot(Trk, Lb(end,:)*Urk-Xr, '.-')