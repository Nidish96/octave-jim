clc
clear all
addpath('../ROUTINES/TRANSIENT/')

%% Parameters
m = 1.0;
k = 1.0;
c = 2.0/(2*sqrt(k*m));

muN0 = 1.0;
kt = 3.0;

function [f, z, dfdx, dfdxd] = eldryfric(t, x, z, xd, muN0, kt)
  f = kt*(x-z);
  dfdx = kt;
  dfdxd = 0;
  if abs(f)>muN0  % Slipped - Update
    f = muN0*sign(f);
    dfdx = 0;
    dfdxd = 0;
    
    z = x - muN0*sign(f)/kt;
  end
end

%% RK Gen
function [Xp, zp] = func(t, X, z, f, m, c, k, muN0, kt)
  [fp, zp] = eldryfric(t, X(1), z, X(2), muN0, kt);
  Xp = [0 1; -k/m -c/m]*X + [0; -fp/m] + [0; f(t)/m];
end

freq = 0.5;
fex = @(t) 5*sin(2*pi*freq*t);

% RK 45 Butcher Tableau
pars.a = [0 0 0 0 0 0; 
          1/4 0 0 0 0 0; 
          3/32 9/32 0 0 0 0; 
          1932/2197 -7200/2197 7296/2197 0 0 0;
          439/216 -8 3680/513 -845/4104 0 0;
         -8/27 2 -3544/2565 1859/4104 -11/40 0];
pars.b = [16/135 0 6656/12825 28561/56430 -9/50 2/55];
pars.c = [0 1/4 3/8 12/13 1 1/2];
% Display
pars.Display = 'on';

t = linspace(0, 40, 1000);
[T, X, z] = RK_GEN_TV(@(t, X, z) func(t, X, z, fex, m, c, k, muN0, kt), t, [0.5; 0], 0, pars);

%% HHTA Hysteretic
ABG = [0, 1/4, 1/2];  % Unconditionally Stable Newmark-Alpha
%% ABG = [0, 1/6, 1/2];  % Implicit linear acceleration
%% ABG = [-0.1, 1/6, 1/2];  % HHT-Alpha

opts = struct('reletol', 1e-6, 'etol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, 'Display', true, 'ITMAX', 10);

[Th, Xh, zh, Xdh, Xddh] = HHTA_NONLIN_HYST(m, c, k, fex,
					   @(t, x, z, xd) eldryfric(t, x, z, xd, muN0, kt),
					   0.5, 0, 0, 0, 40, t(2)-t(1), ABG(1), ABG(2), ABG(3), opts);

%% HBM
Nt = 256;
h = [0 1 2 3 4 5];  Nhc = sum(h==0)+2*sum(h~=0);
Fl = [0 0 5 zeros(1, Nhc-3)]';
U0 = zeros(Nhc,1);

opt = optimset('Jacobian', 'on', 'Display', 'iter');
[U, ~, info, op] = fsolve(@(X) SDOF_NLHYST_HBRESFUN([X; 2*pi*freq], m, c, k, Fl,
						    @(t,x,z,xd) eldryfric(t,x,z,xd,muN0,kt), h, Nt),
			  U0, opt);
th = linspace(0, 2*pi, Nt+1)'/(2*pi*freq); th(end) = [];
ut = TIMESERIES_DERIV(Nt, h, U, 0);

Ncyc = ceil(40*freq);
Thb = reshape(th + [0:Ncyc-1]/freq, [], 1);
uthb = repmat(ut, Ncyc, 1);

figure(1)
clf()
plot(T, X(:, 1), 'k-', T, z, 'k--'); hold on
plot(Th, Xh, 'bx', Th, zh, 'r:');

figure(2)
clf()
plot(T, X(:,1)); hold on
plot(Th, Xh, 'x')
plot(Thb, uthb, 'k--')
