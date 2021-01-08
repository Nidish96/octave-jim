clc
clear all
addpath('../../ROUTINES/')
addpath('../../ROUTINES/HARMONIC')
addpath('../../ROUTINES/SOLVERS')

analyze = false

%% Parameters
m = 1;
k = 1;
c = 2*1e-2;
g = 20;
kg = 3;
n = 2;

h = [0 1 2 3 4 5 6];
Nt = 256;
Nhc = sum(h==0)+2*sum(h~=0);
Fl = zeros(Nhc,1);
Fl(2) = 1;

Fs = [0.01 0.1 0.5 1.2];
Angs = [1e-4 1e-2 1e-1 1e-1];
Wstart = 0.01;
Wend = 2;

Uss = cell(size(Fs));
Ws = cell(size(Fs));
As = cell(size(Fs));
Phs = cell(size(Fs));

if analyze
  for i=1:length(Fs)
    U0 = HARMONICSTIFFNESS(m, c, k, Wstart, h)\Fl;

    ds = 0.01;
    Copt = struct('Nmax', 1000, 'Display', 1, 'angopt', Angs(i), 'dsmax', 0.5, 'dsmin', 0.01);
    Uss{i} = CONTINUE(@(Uw) UNIL_RESFUN(Uw, m, c, k, g, kg, Fl*Fs(i), Nt, h), U0, Wstart, Wend, ds, Copt);

    Ws{i} = Uss{i}(end,:);
%      As{i} = sqrt(sum([1;0.5*ones(Nhc-1,1)].*Uss{i}(1:end-1,:).^2));
    As{i} = sqrt(sum([1;ones(Nhc-1,1)].*Uss{i}(1:end-1,:).^2));    
    Phs{i} = angle((Uss{i}(2, :)-1j*Uss{i}(3, :))./(Fl(2)-1j*Fl(3)));
  end
  save('./DATS/Unildat.mat', 'Uss', 'Fs', 'Ws', 'As', 'Phs', '-7');
else
  load('./DATS/Unildat.mat')
end

%% RQNM
rqnm.x = logspace(-3, 4, 1000);
rqnm.q = rqnm.x*sqrt(m);
rqnm.wp = sqrt((k/m+kg./(rqnm.q*sqrt(m)).*max(rqnm.x-g,0)));
rqnm.wm = sqrt(k/m)*ones(size(rqnm.x));
rqnm.zp = 0.5*c/m./rqnm.wp;
rqnm.zm = 0.5*c/m./rqnm.wm;

%% RQNM2
Nt = 128;
t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
qt = cos(t).*rqnm.q(:)';
xt = qt/sqrt(m);

Lt = k/m+kg./(rqnm.q.*sqrt(m)).*max(rqnm.x-g,0);
Lams = sqrt(sum((GETFOURIERCOEFF(1, Lt.*qt)./rqnm.q(:)').^2))';
xdot = -sin(t).*(rqnm.q(:).*sqrt(Lams))';

fm = (Lt.*xt)/m;
z = GETFOURIERCOEFF(0, xdot.*(Lt.*xt))'./(rqnm.q(:).^2.*Lams.^1.5);

rqnm2.x = rqnm.x;
rqnm2.w = sqrt(Lams);
rqnm2.z = z;

%% RQNM3
Nt = 128;
##opt = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
opt = optimset('Jacobian', 'on', 'Display', 'iter');
x0 = [0; k/m];
sols = zeros(3, length(rqnm.q));
for i=1:length(rqnm.q)
    Q = rqnm.q(i);
    xs = fsolve(@(qw2) RQHBRESFUN([qw2; Q], @(t, q) (k*q/sqrt(m)+max(kg.*(q/sqrt(m)-g),0))/m, ...
        @(t, q) (k/sqrt(m)+(q/sqrt(m)>g).*kg)/m, Nt), x0, opt);
    sols(:, i) = [xs; Q];
    x0 = xs;
    
    fprintf('%d/%d\n', i, length(rqnm.q));
end

rqnm3.q = sols(3, :)';
rqnm3.x = rqnm3.q/sqrt(m);
rqnm3.w = sqrt(sols(2, :));
rqnm3.z = rqnm2.z;

%% Plotting
lw = 2;
fs = 20;
figure(10)
clf()
plot(rqnm.x, rqnm.wp.^2, 'b-', 'LineWidth', lw); hold on
plot(-rqnm.x, rqnm.wm.^2, 'r-', 'LineWidth', lw); hold on
xlim(100*[-1 1])
ylim([0 4])
set(gca, 'fontsize', fs)
xlabel('Response Amplitude (x)')
ylabel('Rayleigh Quotient')
legend('$\lambda^+$', '$\lambda^-$', 'Location', 'northwest')

% print('../PRESENTATION/FIGS/7a_UNILUNSYM.eps', '-depsc')
%%
figure(1)
clf()
for i=1:length(Ws)
  plot(Ws{i}, As{i}, '-', 'LineWidth', lw); hold on
end
aa1 = plot(rqnm.wp, rqnm.x, 'b-', 'LineWidth', lw);
aa2 = plot(rqnm.wm, rqnm.x, 'r-', 'LineWidth', lw);
aa3 = plot(rqnm2.w, rqnm2.x, 'k-.', 'LineWidth', lw*1.5);
aa4 = plot(rqnm3.w, rqnm3.x, 'k-.', 'LineWidth', lw*1.5);
aa4.Visible = false;

xlim([0.5 2])
ylim([0 50])
set(gca, 'fontsize', fs)
xlabel('Forcing Frequency (\omega)', 'Fontsize', fs)
ylabel('Response Amplitude', 'Fontsize', fs)

legend([aa1, aa2, aa3], '$\sqrt{\lambda^+}$', '$\sqrt{\lambda^-}$', 'RQNMA', 'Location', 'southeast')

print('../PRESENTATION/FIGS/3a_UNILRQNM.eps', '-depsc')

aa3.Visible = false;
aa4.Visible = true;
legend([aa1, aa2, aa4], '$\sqrt{\lambda^+}$', '$\sqrt{\lambda^-}$', 'RQNMA', 'Location', 'southeast')

print('../PRESENTATION/FIGS/3b_UNILRQNM.eps', '-depsc')

% plot((rqnm.wp+rqnm.wm)/2, rqnm.x, 'k--', 'LineWidth', lw)
% plot(sqrt(rqnm.wp.*rqnm.wm), rqnm.x, 'k:', 'LineWidth', lw)
% plot(sqrt((rqnm.wp.^2+rqnm.wm.^2)/2), rqnm.x, 'k-.', 'LineWidth', lw)

% print('../PRESENTATION/FIGS/7a_UNILRQNM.eps', '-depsc')
%%
figure(2)
clf()
loglog(rqnm.q, rqnm.zp, 'b-', 'LineWidth', lw); hold on
loglog(rqnm.q, rqnm.zm, 'r-', 'LineWidth', lw);
loglog(rqnm.q, 0.5*c/m./((rqnm.wp+rqnm.wm)/2), 'k--', 'LineWidth', lw);
loglog(rqnm.q, 0.5*c/m./sqrt(rqnm.wp.*rqnm.wm), 'k:', 'LineWidth', lw);
loglog(rqnm.q, 0.5*c/m./sqrt((rqnm.wp.^2+rqnm.wm.^2)/2), 'k-.', 'LineWidth', lw);

xlabel('Modal Amplitude', 'fontsize', fs)
ylabel('Modal Damping Factor', 'fontsize', fs)
set(gca, 'fontsize', fs)
ylim([1e-3 1e-1])

legend('Mode +', 'Mode -', 'Arith. Mean', 'Geometric Mean', 'RMS', 'Location', 'southwest')

% print('../PRESENTATION/FIGS/7b_UNILRQNM_Z.eps', '-depsc')
