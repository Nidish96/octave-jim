clc
clear all

%% Parameters
Ne = 2;
rhoAbL = 7800*(1e-2^2)/(1/Ne);
AEbL = (1e-2^2)*2e11/(1/Ne);

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

fsamp = 2^15;
dt = 1/fsamp;
T0 = 0;
T1 = 0.1;

Ts = T0:dt:T1;

type = 'IMP';
fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw));

% type = 'WGN';
% fext = wgn(size(Ts,1), size(Ts,2), 40);
% fex = @(t) interp1(Ts, fext, t);

sys = @(t, y) [zeros(Ne) eye(Ne);-Mb\Kb -Mb\Cb]*y + [zeros(Ne,1); Mb\Fb]*fex(t) + [zeros(Ne,1); Mb\Fs]*sa;

y0 = [Ustat; Ustat*0];
% tic
% [T, Y] = ode45(sys, Ts, y0); 
% toc

disp('done')

%% Test with MDOFGEN
MDL = MDOFGEN(Mb, Kb, Cb, Lb);
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

% yyaxis right
% plot(T, fex(T), '.-')