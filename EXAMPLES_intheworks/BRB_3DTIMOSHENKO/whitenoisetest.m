clc
clear all
addpath('../../ROUTINES/HARMONIC/')

%% Parameters
Ne = 2;
rhoAbL = 7800*(1e-2^2)/(1/Ne);
AEbL = (1e-2^2)*2e11/(1/Ne);
b = 2e4*0;

Me = [2 1;1 2]*rhoAbL;
Ke = [1 -1;-1 1]*AEbL;

M = zeros(Ne+1);
K = zeros(Ne+1);
for e=1:Ne
    M(e:e+1, e:e+1) = M(e:e+1, e:e+1) + Me;
    K(e:e+1, e:e+1) = K(e:e+1, e:e+1) + Ke;
end
Lb = eye(Ne+1);  Lb(:, 1) = [];

Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;
Fb = Lb(end, :)';

Fs = Lb(end, :)';  % "static force"
sa = 1;
%%
[Vs, Ws] = eig(Kb, Mb);
Ws = sqrt(diag(Ws));

zts = [0.1; 0.2]*1e-2;
ab = [1./(2*Ws(1:length(zts))) Ws(1:length(zts))/2]\zts;

Cb = ab(1)*Mb + ab(2)*Kb;

%% Static Solution
Ustat = Kb\(sa*Fs);

%%
bw = 1000;
famp = 100;

fsamp = 2^16;
dt = 1/fsamp;
T0 = 0;
T1 = 1;

Ts = T0:dt:T1;

type = 'IMP';
fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw));

% type = 'WGN';
% fext = wgn(size(Ts,1), size(Ts,2), 180);
% fex = @(t) interp1(Ts, fext, t);

sys = @(t, y) [zeros(Ne) eye(Ne);-Mb\Kb -Mb\Cb]*y+[zeros(Ne,1); Mb\Fb]*(fex(t)+sa)+[zeros(Ne,1); Mb\(-b*Lb(end,:)'*(Lb(end,:)*y(1:Ne)).^3)];

y0 = [Ustat; Ustat*0];
tic
[T, Y] = ode45(sys, Ts, y0); 
toc
sol = ode45(sys, Ts, y0);
disp('done')

%% Test with MDOFGEN
MDL = MDOFGEN(Mb, Kb, Cb, Lb);
fnl = @(t,u,ud) deal(b*u.^3, 3*b*u.^2, ud*0);
MDL = MDL.SETNLFUN(1+3, Lb(end,:), fnl);

FEX = @(t) fex(t).*Fb + Fs*sa;

opts = struct('Display', 'waitbar');

tic
[TT, U, Ud, Udd, MDL] = MDL.HHTAMARCH(T0, T1, dt, Ustat, zeros(size(Ustat)), ...
    FEX, opts);
toc

%% Test wiith lsim
sssys = ss([zeros(Ne) eye(Ne);-Mb\Kb -Mb\Cb], [zeros(Ne,1); Mb\Fb], kron(eye(2), Lb(end, :)), 0);
tic
Yl = lsim(sssys, fex(Ts)+sa, Ts, y0);
toc

%% Plotting
figure(1)
clf()

% plot(T, Lb(end,:)*Y(:, 1:Ne)', '-'); hold on
plot(Ts, Yl(:,1), '-'); hold on
plot(TT, Lb(end,:)*U, '.'); hold on
xlabel('Time (s)')
ylabel('Response (m)')

figure(2)
clf()

% plot(T, Lb(end,:)*Y(:, Ne+(1:Ne))', '-'); hold on
plot(Ts, Yl(:,2), '-'); hold on
plot(TT, Lb(end,:)*Ud, '.')
xlabel('Time (s)')
ylabel('Velocity (m)')

figure(10)
clf()
plot(T, fex(T), '.-'); hold on
plot(TT, fex(TT), 'o-')
xlim([0 1e-3])

xlabel('Time (s)')
ylabel('Force (N)')

%% FREQ DOM
[fr1, uf1] = FFTFUN(T, (Lb(end,:)*Y(:, 1:Ne)')');
[~, ff1] = FFTFUN(T, fex(T));

[fr2, uf2] = FFTFUN(TT', (Lb(end,:)*U)');
[~, ff2] = FFTFUN(TT, fex(TT'));

[fr3, uf3] = FFTFUN(Ts', Yl(:,1));
[~, ff3] = FFTFUN(Ts', fex(Ts'));

figure(3)
clf()
subplot(2,1,1)
loglog(fr1, abs(uf1./ff1), '-'); hold on
loglog(fr2, abs(uf2./ff2), '.')
loglog(fr3, abs(uf3./ff3), '.')

xlabel('Frequency (Hz)')
ylabel('Response (m)')

subplot(2,1,2)
semilogx(fr1, rad2deg(angle(uf1./ff1)), '-'); hold on
semilogx(fr2, rad2deg(angle(uf2./ff2)), '.')
semilogx(fr3, rad2deg(angle(uf3./ff3)), '.')

xlabel('Frequency (Hz)')
ylabel('Phase (deg)')