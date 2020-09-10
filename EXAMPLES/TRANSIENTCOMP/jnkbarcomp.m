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

% kc = 1e6;
% fnl = @(t,u,ud) deal(kc*u.^3, 3*kc*u.^2, zeros(size(u)));
kt  = 7.5e8;
muN = 7.5e0;
fnl = @(t, u, varargin) JENKFORCE(t, u, kt, muN, varargin{:});
MDL = MDL.SETNLFUN(2+3, Lb(end,:), fnl);

%% March
fsamp = 2^16;
T0 = 0;  T1 = 2000/fsamp;  dt = 1./fsamp;

%% Excitation
bw = 1000;

% IMPULSE
% bw = 2000;
% famp = 1e0;
% type = 'IMP';
% fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw))*famp;
% FEX = @(t) Lb(end,:)'*fex(t);

% White Gaussian Noise
famp = 0.1;
fext = wgn(length(T0:dt:T1), 1, 40+20*log10(famp));
fex = @(t) interp1(T0:dt:T1, fext, t);
FEX = @(t) Lb(end,:)'*fex(t);
%% Solve
opts = struct('Display', 'off');

tic
[Trk, Urk, Udrk, Uddrk] = MDL.RKGENMARCH(T0, T1, dt, zeros(MDL.Ndofs,1), ...
    zeros(MDL.Ndofs,1), FEX, opts);
toc

tic 
[Thh, Uhh, Udhh, Uddhh] = MDL.HHTAMARCH(T0, T1, dt, zeros(MDL.Ndofs,1), ...
    zeros(MDL.Ndofs,1), FEX, opts);
toc

% %% Rate form model [y; ydot; z]
% % func = @(t, y) [zeros(Ne), eye(Ne), zeros(Ne,1); -Mb\Kb, -Mb\Cb, -Mb\MDL.NLTs.L'; zeros(1, 2*Ne+1)]*y + ...
% %     [zeros(2*Ne,1); (abs(y(end))<muN)*(kt*MDL.NLTs.L*y(Ne+(1:Ne)))] + ...
% %     [zeros(MDL.Ndofs,1); Mb\Lb(end,:)'; 0]*fex(t);
% 
% func_stk = @(t, y, x0) [zeros(Ne), eye(Ne); -Mb\Kb, -Mb\Cb]*y + ...
%     [zeros(Ne,1); -Mb\MDL.NLTs.L'*kt*(MDL.NLTs.L*y(1:Ne)-x0)] + ...
%     [zeros(MDL.Ndofs,1); Mb\Lb(end,:)']*fex(t);
% 
% func_slp = @(t, y, x0) [zeros(Ne), eye(Ne); -Mb\Kb, -Mb\Cb]*y + ...
%     [zeros(Ne,1); -Mb\MDL.NLTs.L'*(muN*sign(MDL.NLTs.L*y(Ne+(1:Ne))))] + ...
%     [zeros(MDL.Ndofs,1); Mb\Lb(end,:)']*fex(t);
% 
% % func = @(t, y) [zeros(MDL.Ndofs), eye(MDL.Ndofs); -Mb\Kb, -Mb\Cb]*y + [zeros(MDL.Ndofs,1); Mb\Lb(end,:)']*fex(t);
% 
% tic
% [Tst, yst] = ode45(@(t, y) func_stk(t, y, 0), T0:dt:T1, zeros(Ne*2,1));
% [Tsl, ysl] = ode45(@(t, y) func_slp(t, y, 0), T0:dt:T1, zeros(Ne*2,1));
% toc
% 
% % tic
% % [To, yo] = ode45(func, T0:dt:T1, zeros(Ne*2,1));
% % toc

%% Plots
figure(1)
clf()

yyaxis left
plot(Trk, Lb(end,:)*Urk, 'r.-'); hold on
plot(Thh, Lb(end,:)*Uhh, 'b.-'); hold on
% plot(To, Lb(end,:)*yo(:, 1:MDL.Ndofs)', 'o')

xlabel('Time (s)')
ylabel('Displacement (m)')

% yyaxis right
% plot(Trk, fex(Trk), 'k--')

% ylabel('Excitation')

figure(2)
clf()
plot(Trk, diag(Udrk'*Mb*Udrk+Urk'*Kb*Urk), 'r.-'); hold on
plot(Thh, diag(Udhh'*Mb*Udhh+Uhh'*Kb*Uhh), 'b.-'); hold on

xlabel('Time (s)')
ylabel('Linear Energy')

% yyaxis right
% plot(Trk, fex(Trk), 'k--')

% ylabel('Excitation')
% grid on

%% Frequency Domain
[freqs, Ff] = FFTFUN(Thh(1:end-1)', fex(Thh(1:end-1))');
[~, Uhf] = FFTFUN(Thh(1:end-1), (Lb(end,:)*Uhh(:, 1:end-1))');
[~, Urf] = FFTFUN(Thh(1:end-1), interp1(Trk, (Lb(end,:)*Urk)', Thh(1:end-1))');

figure(3)
% clf()
% plot(freqs, abs(Urf./Ff), '.-'); hold on
plot(freqs, abs(Uhf./Ff), '.-'); hold on

xlabel('Frequency (Hz)')
ylabel('Response Amplitude')

set(gca, 'YScale', 'log')