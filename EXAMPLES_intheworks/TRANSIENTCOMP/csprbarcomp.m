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

% Nonlinearity
bt = 1e4;
cspr = @(t, u, ud) deal(bt*u.^3, 3*bt*u.^2, zeros(size(u)));
MDL = MDL.SETNLFUN(1+3, Lb(end,:), cspr);

%% Excitation
bw = 1000;

% IMPULSE
bw = 1000;
famp = 1e12;
type = 'IMP';
fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw));
FEX = @(t) Lb(end,:)'*fex(t)*famp;

%% March
fsamp = 2^16;
T0 = 0;  T1 = 400/fsamp;  dt = 1./fsamp;

opts = struct('Display', 'waitbar');
tic
[Trk, Urk, Udrk, Uddrk] = MDL.RKGENMARCH(T0, T1, dt, zeros(MDL.Ndofs,1), ...
    zeros(MDL.Ndofs,1), FEX, opts);
toc

tic
[Thh, Uhh, Udhh, Uddhh] = MDL.HHTAMARCH(T0, T1, dt, zeros(MDL.Ndofs,1), ...
    zeros(MDL.Ndofs,1), FEX, opts);
toc
%% 
func = @(t, y) [zeros(Ne),eye(Ne); -Mb\Kb,-Mb\Cb]*y + [zeros(Ne,1); Mb\Lb(end,:)']*fex(t)*famp + [zeros(Ne,1); -Mb\MDL.NLFORCE(t, y(1:Ne), y(Ne+(1:Ne)), 0)];
tic
[Tho, yho] = ode45(func, Thh, zeros(Ne*2,1));
toc

tic
[Tro, yro] = ode45(func, Trk, zeros(Ne*2,1));
toc
%% Plots
figure(1)
clf()

plot(Trk, Lb(end,:)*Urk, '.-'); hold on
plot(Thh, Lb(end,:)*Uhh, '+-'); hold on
plot(Tho, Lb(end,:)*yho(:, 1:Ne)', 'o')

xlabel('Time (s)')
ylabel('Displacement (m)')

yyaxis right
plot(Trk, fex(Trk), 'k--')

ylabel('Excitation')

figure(2)
clf()
plot(Thh, Lb(end,:)*(yho(:,1:Ne)'-Uhh), 'b.-'); hold on
plot(Trk, Lb(end,:)*(yho(:,1:Ne)'-Urk), 'r.-'); hold on
% plot(To, interp1(Thh, Lb(end,:)*Uhh, To)-Lb(end,:)*yo(:, 1:Ne)', 'b.')
% plot(To, interp1(Trk, Lb(end,:)*Urk, To)-Lb(end,:)*yo(:, 1:Ne)', 'r.-')