clc
clear all

addpath('../../ROUTINES/')
addpath('../../ROUTINES/SOLVERS/')
addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/FEM/')
addpath('../../ROUTINES/FEM/BEAMS/')

Nein = 8;  % Number of Elements in the Bolt Area
load(sprintf('./DATA/%dIN_MATS.mat',Nein), 'M', 'K', 'Fbolt', 'L1', ...
     'L2', 'SensorLocs', 'RECOV', 'R1', 'R2', 'BM1', 'IN1', 'Beam1', ...
     'BM2', 'IN2', 'Beam2', 'pars', 'parsint', 'Nebb', 'Nein', ...
     'wdt', 'nu', 'int1nds', 'int2nds', 'remnds', ...
     'i1s', 'i2s', 'Vrbms');

Lrbms = null(full(Vrbms'*M));
Lrbms(abs(Lrbms)<1e-10) = 0;
Lrbms = sparse(Lrbms);


LML = Lrbms'*M*Lrbms;
LKL = Lrbms'*K*Lrbms;
LFb = Lrbms'*Fbolt;

% Truncate sparse
LML(abs(LML)<1e-10) = 0;
LKL(abs(LKL)<1e-10) = 0;

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
LTrel(abs(LTrel)<1e-10) = 0;
QrelL(abs(QrelL)<1e-10) = 0;

LTrel = sparse(LTrel);
QrelL = sparse(QrelL);

%% Contact Model
Aint = sum(sum(T1(1:6:end, :)));

Pint = Prestress/Aint;
sint = 100e-6;
chi  = 2.0;
ktkn = chi*(1-nu)/(2-nu);
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn   = kt/ktkn;

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

fnl = @(t, u, varargin) ELDRYFRICT3D(t, u, kt, kt, kn, mu, gap, varargin{:});
% fnl = @(t, u, varargin) LIN3D(t, u, kt, kt, kn, mu, gap, varargin{:});
MDL = MDL.SETNLFUN(2+5, QrelL, fnl, LTrel);

%% Static Prestress Analysis
U0 = (LKL+Ks)\(LFb*Prestress);
% [R, J] = MDL.STATRESFUN(U0, LFb*Prestress);

opts = struct('reletol', 1e-6, 'Display', true);
[Ustat, ~, ~, ~, J0] = NSOLVE(@(U) MDL.STATRESFUN(U, LFb*Prestress), U0, opts);

Ws = sort(sqrt(eigs(J0, LML, 10, 'SM')));
disp(Ws/2/pi)

% fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% Ustat = fsolve(@(U) MDL.STATRESFUN(U, LFb*Prestress), U0*10, fopts);

% disp(sqrt(eigs(J0, LML, 10, 'SM'))/2/pi);

%% HBM Fresp
h = 0:1;
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^7;
Nd = size(LML,1);

fa = 10.0;

if Nhc>1
    Fl = kron([0; 1; 0; zeros(Nhc-3,1)], fa*Lrbms'*RECOV(3,:)');
else
    Fl = zeros(Nd, 1);
end
Fl(1:Nd) = LFb*Prestress;  % Static Solution

Wst = 100*2*pi;
Wen = 200*2*pi;

Wst = 200*2*pi;
Wen = 100*2*pi;
dw = 1*2*pi;

E = HARMONICSTIFFNESS(MDL.M, MDL.C, MDL.K+J0, Wst, h);
U0 = E\Fl;
U0(1:Nd) = Ustat;

% opts = struct('reletol', 1e-6, 'etol', 1e-12, 'utol', 1e-12, 'Display', true);
% Us = NSOLVE(@(U) MDL.HBRESFUN([U; Wst], Fl, h, Nt, 1e-6), U0, opts);

Copt = struct('Nmax', 1000, 'itDisplay', false);
% Copt.dsmax = dw;

UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

figure(1)
% clf()

plot(UC(end,:)/2/pi, sqrt([1 0.5*ones(1, Nhc-1)]*(kron(eye(Nhc), RECOV(3,:)*Lrbms)*UC(1:end-1,:)).^2)/fa, '.-'); hold on
plot(Ws(1)*[1 1]/2/pi, ylim, 'k--');
xlabel('Frequency (Hz)')
ylabel('RMS Response (m)')

%% NLvib continuation
Sopt = struct('jac', 'full', 'stepmax', 100, 'MaxFfunEvals', 100);

Uc = solve_and_continue(U0, @(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), Wst, Wen, dw, Sopt);
