% clc
% clear all

addpath('../../ROUTINES/')
addpath('../../ROUTINES/SOLVERS/')
addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/FEM/')
addpath('../../ROUTINES/FEM/BEAMS/')

Nein = 8;  % Number of Elements in the Bolt Area
load(sprintf('./MATS/%dIN_MATS.mat',Nein), 'M', 'K', 'Fbolt', 'L1', ...
     'L2', 'SensorLocs', 'RECOV', 'R1', 'R2', 'BM1', 'IN1', 'Beam1', ...
     'BM2', 'IN2', 'Beam2', 'pars', 'parsint', 'Nebb', 'Nein', ...
     'wdt', 'nu', 'int1nds', 'int2nds', 'remnds', ...
     'i1s', 'i2s', 'Vrbms');

Lrbms = null(full(Vrbms'*M));
Lrbms(abs(Lrbms)<eps) = 0;
Lrbms = sparse(Lrbms);


LML = Lrbms'*M*Lrbms;
LKL = Lrbms'*K*Lrbms;
LFb = Lrbms'*Fbolt;

% Truncate sparse
LML(abs(LML)<eps) = 0;
LKL(abs(LKL)<eps) = 0;

Prestress = 12e3;
%% Set up Quadrature
No = 2;

Les = diff(IN1.X);
Wys = IN1.WY(1:end-1);
Zs  = 0-IN1.Z(1:end-1);
[Q1, T1] = TM3D_ND2QP(Les, Wys, Zs, No);

Les = diff(IN2.X);
Wys = IN2.WY(1:end-1);
Zs  = 0-IN2.Z(1:end-1);
[Q2, T2] = TM3D_ND2QP(Les, Wys, Zs, No);

Qrel = zeros(Nein*No^2*3, (Nebb+Nein+1)*2*6);
Trel = zeros((Nebb+Nein+1)*2*6, Nein*No^2*3);

Qrel(:, i1s) = Q1;
Qrel(:, i2s) = -Q2;

Trel(i1s, :) = T1;
Trel(i2s, :) = -T2;

LTrel = Lrbms'*Trel;
QrelL = Qrel*Lrbms;

% Truncate for Storage Efficiency
% LTrel(abs(LTrel)<1e-10) = 0;
% QrelL(abs(QrelL)<1e-10) = 0;

LTrel = sparse(LTrel);
QrelL = sparse(QrelL);

%% Contact Model
Aint = sum(sum(T1(1:6:end, :)));

Pint = Prestress/Aint;
sint = 1e-6;
chi  = 2.0;
ktkn = chi*(1-nu)/(2-nu);
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn   = kt/ktkn

kn = 1e10;
% kt = kn*ktkn;

%% Linear Dissipation
tstiff = kron(ones(Nein*No^2,1), [kt;kt;kn]);
Ks = LTrel*diag(tstiff)*QrelL;

[Vs, Ws] = eigs(LKL+Ks, LML, 10, 'SM');
[Ws, si] = sort(sqrt(diag(Ws)));
Vs = Vs(:, si);
Vs = Vs./sqrt(diag(Vs'*LML*Vs))';

zetas = [0.1; 0.2]*1e-2;  % 0.1%, 0.2%
ab = [1./(2*Ws(1:length(zetas))) Ws(1:length(zetas))/2]\zetas;

LCL = ab(1)*LML + ab(2)*(LKL+Ks);

%% Setup Class Object
MDL = MDOFGEN(LML, LKL, LCL, Lrbms);

mu   = 0.25;
% mu = 1e6;
gap  = 0;

% fnl = @(t, u, varargin) ELDRYFRICT3D(t, u, kt, kt, kn, mu, gap, varargin{:});
fnl = @(t, u, varargin) LIN3D(t, u, kt, kt, kn, mu, gap, varargin{:});
MDL = MDL.SETNLFUN(2+5, QrelL, fnl, LTrel);
% MDL = MDL.SETNLFUN(2+5, QrelL, fnl, QrelL');

%% Static Prestress Analysis
U0 = (LKL+Ks)\(LFb*Prestress);
% [R, J] = MDL.STATRESFUN(U0, LFb*Prestress);

opts = struct('reletol', 1e-6, 'Display', true);
[Ustat, ~, ~, ~, J0] = NSOLVE(@(U) MDL.STATRESFUN(U, LFb*Prestress), U0, opts);

[V, Ws] = eig(full(J0), full(LML));
[Ws, si] = sort(sqrt(diag(Ws)));
V = V(:, si);
disp(max(Ws/2/pi))

%% Excitation
fsamp = 2^18;  % Sampling frequency (2^18)
T0 = 0;  T1 = 400/fsamp;  dt = 1/fsamp;

%% ldof = 6;
%% DOF = 'Z'
% ldof = 8;
% DOF = 'Y'
ldof = 1;
DOF = 'X'

% % IMPULSE
% bw = 1000;
% famp = 100;
% type = 'IMP';
% fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw))*famp;
% fext = fex(T0:dt:T1);
% FEX = @(t) Lrbms'*RECOV(ldof,:)'*fex(t)+LFb*Prestress;

% WHITE GAUSSIAN NOISE
rng(1)
bw = -1;
famp = 1
type = 'WGN';
fex = @(t) wgn(size(t,1), size(t,2), 40+20*log10(famp));
fext = fex(T0:dt:T1);
FEX = @(t) Lrbms'*RECOV(ldof,:)'*interp1(T0:dt:T1, fext, t)+LFb*Prestress;

[freqs, Ff] = FFTFUN((0:(1/fsamp):1)', fex(0:(1/fsamp):1)');

% figure(10)
% clf()
% semilogy(freqs, abs(Ff), '.'); hold on
% for i=1:length(Ws)
%     semilogy(Ws(i)*[1 1]/2/pi, ylim, 'k--')
% end
% xlim([0 3e3])
% xlabel('Frequency (Hz)')
% ylabel('Forcing (N)')

%% March
opts = struct('Display', 'waitbar', 'minstep', 1e-8);

% [~, ~, ~, MDL] = MDL.NLFORCE(0, Ustat, zeros(size(Ustat)), 0, 1);

[Trk, Urk, Udrk, Uddrk, ~] = MDL.RKGENMARCH(T0, T1, dt, Ustat, ...
    zeros(size(Ustat)), FEX, opts);

opts = struct('Display', 'waitbar');
[Thh, Uhh, Udhh, Uddhh, ~] = MDL.HHTAMARCH(T0, T1, dt, Ustat, ...
    zeros(size(Ustat)), FEX, opts);

%% Linearized System (analytical)
[~, JNL] = MDL.NLFORCE(0, Ustat, zeros(size(Ustat)), 0);
[~, Tha, Xha] = lsim(ss([zeros(MDL.Ndofs) eye(MDL.Ndofs); -MDL.M\(MDL.K+JNL) -MDL.M\MDL.C], ...
                [zeros(MDL.Ndofs,2); MDL.M\[Lrbms'*RECOV(ldof,:)' LFb*Prestress]], [RECOV(ldof, :)*Lrbms zeros(1, MDL.Ndofs)], ...
                [0 0]), [interp1(T0:dt:T1, fext, Thh); ones(size(Thh))]', Thh, [Ustat; Ustat*0]);

% %% 
% MDlin = MDOFGEN(MDL.M, MDL.K+JNL, MDL.C, MDL.L);
% [Thl, Uhl, Udhl, Uddhl, ~] = MDlin.HHTAMARCH(T0, T1, dt, Ustat, ...
%     zeros(size(Ustat)), FEX, opts);

%% 
odof = 3;

figure(1)
clf()

plot(Trk, RECOV(odof, :)*Lrbms*Urk, 'b.-'); hold on
plot(Thh, RECOV(odof, :)*Lrbms*Uhh, 'ro-'); hold on
% plot(Thl, RECOV(odof, :)*Lrbms*Uhl, 'cx-'); hold on
plot(Thh, RECOV(odof, :)*Lrbms*Xha(:, 1:MDL.Ndofs)', 'k--'); hold on

xlabel('Time (s)')
ylabel('Output disp (m)')
% yyaxis right

% plot(Tha, RECOV(odof, :)*Lrbms*Xha(:, 1:MDL.Ndofs)', 'k--'); hold on

% plot(Thh, fex(Thh), 'k--')












