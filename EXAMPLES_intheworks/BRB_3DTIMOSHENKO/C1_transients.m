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

% Lrbms = null(full(Vrbms'*M));
% % Lrbms(abs(Lrbms)<1e-10) = 0;
% % Lrbms = sparse(Lrbms);
% 
% 
% LML = Lrbms'*M*Lrbms;
% LKL = Lrbms'*K*Lrbms;
% LFb = Lrbms'*Fbolt;

% Truncate sparse
% LML(abs(LML)<1e-10) = 0;
% LKL(abs(LKL)<1e-10) = 0;

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

% LTrel = Lrbms'*Trel;
% QrelL = Qrel*Lrbms;

% Truncate for Storage Efficiency
% LTrel(abs(LTrel)<1e-10) = 0;
% QrelL(abs(QrelL)<1e-10) = 0;

% LTrel = sparse(LTrel);
% QrelL = sparse(QrelL);

%% Contact Model
Aint = sum(sum(T1(1:6:end, :)));

Pint = Prestress/Aint;
sint = 1e-6;
chi  = 2.0;
ktkn = chi*(1-nu)/(2-nu);
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn   = kt/ktkn;
kn = 1e10;

%% Linear Dissipation
tstiff = kron(ones(Nein*No^2,1), [kt;kt;kn]);
Ks = Trel*diag(tstiff)*Qrel;

[Vs, Ws] = eigs(K+Ks, M, 20, 'SM');
[Ws, si] = sort(sqrt(diag(Ws)));
Vs = Vs(:, si);
Vs = Vs./sqrt(diag(Vs'*M*Vs))';

zetas = [0.2; 0.1; 0.4]*1e-2;  % 0.1%, 0.2%
ab = [1./(2*Ws(6+(1:length(zetas)))) Ws(6+(1:length(zetas)))/2]\zetas;

C = ab(1)*M + ab(2)*(K+Ks);

diag(Vs(:,7:end)'*C*Vs(:,7:end))./(2*Ws(7:end))

%% Null space of RBMs
Vrbms = Vrbms(:, 1:6);
Lrbms = null(full(Vrbms'*M));

LML = Lrbms'*M*Lrbms;
LKL = Lrbms'*K*Lrbms;
LCL = Lrbms'*C*Lrbms;
LFb = Lrbms'*Fbolt;

QrelL = Qrel*Lrbms;
LTrel = Lrbms'*Trel;

Ks = Lrbms'*Ks*Lrbms;

%% Setup Class Object
MDL = MDOFGEN(LML, LKL, LCL, Lrbms);

mu   = 0.25;
% mu = 1e6;
gap  = 0;

fnl = @(t, u, varargin) ELDRYFRICT3D(t, u, kt, kt, kn, mu, gap, varargin{:});
% fnl = @(t, u, varargin) LIN3D(t, u, kt, kt, kn, mu, gap, varargin{:});
MDL = MDL.SETNLFUN(2+5, QrelL, fnl, LTrel);
% MDL = MDL.SETNLFUN(2+5, QrelL, fnl, QrelL');

%% Static Prestress Analysis
U0 = (LKL+Ks)\(LFb*Prestress);
% [R, J] = MDL.STATRESFUN(U0, LFb*Prestress);

opts = struct('reletol', 1e-6, 'Display', true);
Ustat = NSOLVE(@(U) MDL.STATRESFUN(U, LFb*Prestress), U0, opts);

[~, J0, MDL] = MDL.STATRESFUN(Ustat, LFb*Prestress);

[V, Ws] = eigs(J0, LML, 10, 'SM');
[Ws, si] = sort(sqrt(diag(Ws)));
V = V(:, si);
disp(Ws/2/pi)

%% Excitation
bw = 1000;
famp = 600;
fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw));
FEX = @(t) famp*Lrbms'*RECOV(3,:)'*fex(t)+LFb*Prestress;
% FEX = @(t) famp*Lrbms'*RECOV(end,:)'*fex(t)+LFb*Prestress;

% FEX = @(t) famp*Lrbms'*sum(RECOV(1:3,:))'*fex(t)+LFb*Prestress;

[freqs, Ff] = FFTFUN((0:1e-4:0.1)', fex(0:1e-4:0.1)');

figure(10)
clf()
semilogy(freqs, abs(Ff)); hold on
for i=1:length(Ws)
    semilogy(Ws(i)*[1 1]/2/pi, ylim, 'k--')
end
xlabel('Frequency (Hz)')
ylabel('Forcing (N)')

%% HHTA
opts = struct('Display', 'waitbar');

T0 = 0;
T1 = 0.1;
dt = 5e-6;

MDL.NLTs.fp = 0*MDL.NLTs.fp;
MDL.NLTs.up = 0*MDL.NLTs.up;

[T, U, Ud, Udd, MDL] = MDL.HHTAMARCH(T0, T1, dt, Ustat, zeros(size(Ustat)), ...
    FEX, opts);

%%
figure(1)
% clf()

% plot(T, RECOV(end-2,:)*Lrbms*U); hold on
% plot(T, RECOV(end-1,:)*Lrbms*U)
plot(T, RECOV(end,:)*Lrbms*Udd/famp); hold on

xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
%%
figure(4)
% clf()
[freqs, Uf] = FFTFUN(T(:), (RECOV(end, :)*Lrbms*Udd)');
[~, Ff] = FFTFUN(T(:), fex(T(:)));
semilogy(freqs, abs(Uf./Ff)/famp)
% semilogy(freqs, abs(Uf))
hold on
for i=1:length(Ws)
    semilogy(Ws(i)*[1 1]/2/pi, ylim, 'k--')
end
xlabel('Frequency (Hz)')
ylabel('FRF Amplitude')
xlim([0 3e3])

% %%
% timei = 0.006;
% 
% % for timei = 0:0.0001:0.01
% ti = 1+fix(timei/dt);
% sc = 1e6/famp;
% 
% sis = [1 size(SensorLocs,1)];
% 
% soln = Lrbms*(U(:, ti)-Ustat);
% % soln = Lrbms*V(:, 1);
% % soln = Vrbms(:, 6);
% 
% sdat = reshape(RECOV*soln, 3, [])';
% 
% figure(5)
% clf()
% 
% DEPICTBEAM_TM3D(diff(Beam1.X), Beam1.WY, Beam1.WZ, ...
%     [Beam1.X, Beam1.Y, Beam1.Z], L1*soln*sc, 'b', ...
%     0.1, 2); hold on
% DEPICTBEAM_TM3D(diff(Beam2.X), Beam2.WY, Beam2.WZ, ...
%     [Beam2.X, Beam2.Y, Beam2.Z], L2*soln*sc, 'r', ...
%     0.1, 2);
% 
% plot3(SensorLocs(sis,1)+sc*sdat(sis, 1), SensorLocs(sis,2)+sc*sdat(sis, 2), ...
%     SensorLocs(sis,3)+sc*sdat(sis, 3), 'k.', 'MarkerFaceColor', 'k', ...
%     'MarkerSize', 40)
% 
% axis equal
% grid on
% 
% xlabel('X Coordinate')
% ylabel('Y Coordinate')
% zlabel('Z Coordinate')
% 
% title(sprintf('Frame %d', ti))
% 
% % pause(0.1)
% % end