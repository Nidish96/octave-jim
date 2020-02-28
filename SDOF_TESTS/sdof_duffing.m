clc
clear all
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/HARMONIC/')

%% Parameters
m = 1.0;
k = 4.0;
c = 2.0/(2*sqrt(k*m));
b = 0.10;

function [f, dfdx, dfdxd] = duffnl(t, x, xd, b)
  f = b*x.^3;
  dfdx = 3*b*x.^3;
  dfdxd = x*0;
end

%% ODE45
func = @(t, X, f) [0 1; -k/m -c/m]*X + [0; -duffnl(t, X(1), X(2), b)/m] + [0; f(t)/m];
freq = 0.5;
fex = @(t) sin(2*pi*freq*t);

t = linspace(0, 40, 1000);
[T, X] = ode45(@(t, X) func(t, X, fex), t, [0.5; 0]);

%% HHTA NONLIN
ABG = [0, 1/4, 1/2];  % Unconditionally Stable Newmark-Alpha
%% ABG = [0, 1/6, 1/2];  % Implicit linear acceleration
%% ABG = [-0.1, 1/6, 1/2];  % HHT-Alpha

opts = struct('reletol', 1e-6, 'etol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, 'Display', true, 'ITMAX', 10);

[Th, Xh, Xdh, Xddh] = HHTA_NONLIN(m, c, k, fex, @(t, x, xd) duffnl(t, x, xd, b),
				  0.5, 0, 0, 40, t(2)-t(1), ABG(1), ABG(2), ABG(3), opts);

%% HBM
Nt = 256;
h = [0 1 2 3 4 5];  Nhc = sum(h==0)+2*sum(h~=0);
Fl = [0 0 1 zeros(1, Nhc-3)]';
U0 = zeros(Nhc,1);

opt = optimset('Jacobian', 'on', 'Display', 'iter');
U = fsolve(@(X) SDOF_NL_HBRESFUN([X; 2*pi*freq], m, c, k, Fl, @(t,x,xd) duffnl(t,x,xd,b), h, Nt),
	   U0, opt);
th = linspace(0, 2*pi, Nt+1)'/(2*pi*freq); th(end) = [];
ut = TIMESERIES_DERIV(Nt, h, U, 0);

Ncyc = ceil(40*freq);
Thb = reshape(th + [0:Ncyc-1]/freq, [], 1);
uthb = repmat(ut, Ncyc, 1);

figure(1)
clf()
plot(T, X(:, 1)); hold on;
plot(Th, Xh, 'x')
plot(Thb, uthb, 'k--')
grid on
