clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/export_fig/')

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if ~isOctave
  set(0,'defaultAxesTickLabelInterpreter', 'default');
  set(0,'defaultTextInterpreter','latex');
  set(0, 'DefaultLegendInterpreter', 'latex');
  set(0,'defaultAxesFontSize',13);
end

analyze = false;
plotout = true;
%% Parameters
M = eye(2);
K = [2 -1;-1 2];
C = 0.01*K;
b = 0.5;

MDL = MDOFGEN(M, K, C, eye(2));
MDL = MDL.SETNLFUN(1+3, [1 0], @(t, u, ud) deal(b*u.^3, 3*b*u.^2, ud*0));

[V, D] = eig(K, M);
[D, si] = sort(diag(D));
V = V(:,si);
%% EPMC
h = [0 1 2 3 4 5 6 7];
% h = 1:7;
Nt = 512;
Nhc = sum(h==0)+2*sum(h~=0);
Fl = zeros(Nhc*2,1);
Fl(sum(h==0)+2) = 1;

Astart = -2;
Aend = 4;
da = 0.01;
% Copt = struct('Nmax', 2000, 'Display', 1, 'DynDscale', 1, 'dsmin', ...
%                 0.0001, 'angopt', 1e-1, 'dsmax', 0.5);
Copt = struct('Nmax', 2000, 'Display', 1, 'DynDscale', 1, 'angopt', 1e-1, ...
    'solverchoice', 2);
% Copt.Dscale = [kron([1e-8; 1e-1*ones(Nhc-1,1)],[1;1]); 1.0; 1e-8; 1.0];
% Copt.Dscale = [ones(Nhc*2,1)*1e-1; 1; 1e-2; 1.0];

BB = cell(2,1);
if analyze
    U0 = kron([zeros(h(1)==0,1); 1; 1; zeros(Nhc-2-sum(h==0),1)], V(:, 1));
    BB{1} = CONTINUE(@(Uwxa) MDL.EPMCRESFUN(Uwxa, Fl, h, Nt, 1e-6), [U0; sqrt(D(1)); 0.01], Astart, Aend, da, Copt);
    
    U0 = kron([0; 1; 1; zeros(Nhc-3,1)], V(:, 2));
    BB{2} = CONTINUE(@(Uwxa) MDL.EPMCRESFUN(Uwxa, Fl, h, Nt, 1e-6), [U0; sqrt(D(2)); 0.01], Astart, Aend, da, Copt);
    
    save('./DATA/2dofdat.mat', 'BB');    
else
    load('./DATA/2dofdat.mat', 'BB');    
end

%% Take only Monotonic Portion
pks = -findpeaks(-BB{1}(end,:));
mi = find((BB{1}(end,1:end-1)-pks(1)).*(BB{1}(end,2:end)-pks(1))<=0, 1);

Psis = (BB{1}(MDL.Ndofs+(1:MDL.Ndofs), [1 1:mi])-1j*BB{1}(2*MDL.Ndofs+(1:MDL.Ndofs), [1 1:mi])); % First Harmonic Mode Shape
qs = [eps 10.^BB{1}(end, 1:mi)];  % Modal Amplitude
w0s = BB{1}(end-2, [1 1:mi]);  % Natural Frequencies
xis = BB{1}(end-1, [1 1:mi]);  % Damping Coefficients

% Psis = Psis(:, 2:end);
% qs = qs(2:end);
% w0s = w0s(2:end);
% xis = xis(2:end);

%% Transient Ringdown
q0 = qs(35);
bt0 = 0;
psi0 = interp1(qs, Psis.', q0).';
w0 = interp1(qs, w0s, q0);

u0 = real(q0*psi0*exp(1j*bt0));
ud0 = -w0*imag(q0*psi0*exp(1j*bt0));

tvec = 0:1e-1:800;

%% Analytical Model
% dqbt = @(t, qbt) 1/(2*interp1(log10(qs), w0s, log10(qbt(1))))*[-interp1(log10(qs), xis, log10(qbt(1)))*interp1(log10(qs), w0s, log10(qbt(1)))*qbt(1); 0];
% [t, qbtt] = ode45(dqbt, tvec, [q0;bt0]);
% [~, phit] = ode45(@(tt, phi) interp1(t, interp1(log10(qs), w0s, log10(qbtt(:,1))), tt), tvec, 0);
% ansol = qbtt(:,1).*interp1(log10(qs), Psis.', log10(qbtt(:,1))).*exp(1j*(phit + qbtt(:,2)));

dqbpht = @(t, qbph) [1/(2*interp1(log10(qs), w0s, log10(qbph(1))))*[-interp1(log10(qs), xis, log10(qbph(1)))*interp1(log10(qs), w0s, log10(qbph(1)))*qbph(1); 0]; interp1(log10(qs), w0s, log10(qbph(1)))];
[t, qbpht] = ode45(dqbpht, tvec, [q0;bt0;0]);

ansol = qbpht(:,1).*interp1(log10(qs), Psis.', log10(qbpht(:,1))).*exp(1j*(qbpht(:,3) + qbpht(:,2)));

ua = real(ansol);
uda = imag(1j*interp1(log10(qs), w0s, log10(qbpht(:,1))).*ansol);

uA = abs(qbpht(:, 1).*interp1(log10(qs), Psis.', log10(qbpht(:,1))));

%% Multiple-Scales Derivation
dqwdt = @(t, qw) -0.5*interp1(log10(qs.*w0s), xis, log10(qw))*qw;

[t, qwt] = ode45(dqwdt, tvec, q0*w0);
wt = interp1(log10(qs.*w0s), w0s, log10(qwt));
qt = interp1(log10(qs.*w0s), qs, log10(qwt));
[~, bt] = ode45(@(tm,bt) interp1(t, wt, tm), tvec, bt0);

phit = interp1(log10(qs.*w0s), Psis.', log10(qwt));  % Mode shapes
ut = real(qt.*phit.*exp(1j*bt));  % solution in physical coordinates

%% Numerical Model
opt = struct('Display', 'waitbar');
[T, U, Ud, Udd] = MDL.HHTAMARCH(0, tvec(end), tvec(2), u0, ud0, @(t) zeros(2,1), opt);

%% Plot Numerical and Analytical Resultsa
figure(2)
clf()
subplot(2,2,1)
plot(T, U(1,:), 'k-'); hold on
plot(t, ua(:,1), 'b-.'); hold on
xlabel('Time (s)')
ylabel('DOF 1')
legend('Reference', 'KSW', 'Location', 'southeast')
title('Krack, Scheidt, Wallaschek (2014) [KSW]')
subplot(2,2,2)
plot(T, U(1,:), 'k-'); hold on
plot(t, ut(:,1), 'r-.')
xlabel('Time (s)')
ylabel('DOF 1')
legend('Reference', 'MS', 'Location', 'southeast')
title('Multiple Scales Formulation (Eq. (26)) [MS]')
subplot(2,2,3)
semilogx(10.^BB{1}(end, :), BB{1}(end-2,:), '.-'); hold on
semilogx(q0, w0, 'ro', 'MarkerFaceColor', 'r')
xlabel('Modal Amplitude')
ylabel('Natural Frequency')
subplot(2,2,4)
plot(T, ua(:,1)'-U(1,:), 'b-'); hold on
plot(T, ut(:,1)'-U(1,:), 'r-'); 
legend('KSW', 'MS', 'Location', 'southeast')
xlabel('Time (s)')
ylabel('Error $(u_{an}-u_{ref})$')
% plot(T, U(2,:), 'b-'); hold on
% plot(t, ua(:,2), 'r-')
% plot(t, uA(:,2), 'k--', 'LineWidth', 2)
% xlabel('Time (s)')
% ylabel('DOF 2')
set(gcf, 'Color', 'white')
export_fig('./FIGS/slowdynobs.png', '-dpng')