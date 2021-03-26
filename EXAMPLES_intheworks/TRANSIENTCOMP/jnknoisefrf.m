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
T0 = 0;  T1 = 2^12/fsamp;  dt = 1./fsamp;

%% Excitation
bw = 1000;

% IMPULSE
% bw = 2000;
% famp = 1e0;
% type = 'IMP';
% fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw))*famp;
% FEX = @(t) Lb(end,:)'*fex(t);

% White Gaussian Noise
famp = 0.05;  % 0.1 
reps = 16;
fext = wgn(length(T0:dt:T1), reps, 40+20*log10(famp));
fex = @(t, r) interp1(T0:dt:T1, fext(:, r), t);
FEX = @(t, r) Lb(end,:)'*fex(t, r);

%% Solve for each repeat
opts = struct('Display', 'off');
ures = zeros(size(fext));
for r=1:reps
    [Thh, Uhh, Udhh, Uddhh] = MDL.HHTAMARCH(T0, T1, dt, zeros(MDL.Ndofs,1), ...
        zeros(MDL.Ndofs,1), @(t) FEX(t, r), opts);
    
    ures(:, r) = Lb(end, :)*Uhh;    
    fprintf('Done %d\n', r);
end

%% FRF estimate
% wndw = ones(size(Thh));  % Rectangular
wndw = hanning(length(Thh), 'periodic');

[freqs, Ffs] = FFTFUN(Thh(:), (fext.*wndw));
[~, Ufs] = FFTFUN(Thh(:), (ures.*wndw));

% % Hanning window in frequency domain directly
% [freqs, Ffs] = FFTFUN(Thh(:), fext);
% [~, Ufs] = FFTFUN(Thh(:), ures);
% 
% freqs = freqs(1:end-1)+0.5*freqs(2);
% Ffs = diff(Ffs);
% Ufs = diff(Ufs);

figure(3)
% clf()

% semilogy(freqs, abs(Ufs./Ffs), '.'); hold on
% semilogy(freqs, mean(abs(Ufs./Ffs), 2), '.'); hold on

semilogy(freqs, abs(mean(Ufs.*conj(Ffs), 2)./mean(Ffs.*conj(Ffs), 2)), '-'); hold on
% semilogy(freqs, abs(mean(Ufs.*conj(Ufs), 2)./mean(Ffs.*conj(Ufs), 2)), 'r-')
legend(sprintf('F=%.2f N', famp))

xlabel('Frequency (Hz)')
ylabel('FRF Amplitude (m/N)')
