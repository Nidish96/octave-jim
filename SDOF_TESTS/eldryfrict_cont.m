clc
clear all
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

%% Parameters
m = 1.0;
k = 1.0;
c = 0.10*2*sqrt(k/m);
% c = 2.0/(2*sqrt(k*m));

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

%% HBM
Nt = 256;
h = [0 1]; Nhc = sum(h==0)+2*sum(h~=0);
Fl = [0 1 0 zeros(1, Nhc-3)]';
U0 = zeros(Nhc,1);
freq = 2.0;

FF = [0.01 0.25 1.50];


opts = struct('reletol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, 'etol', ...
              1e-6, 'ITMAX', 100, 'Display', true);
Copt = struct('Nmax', 100, 'Display', 1, 'angopt', 1e-2, 'opt', opts);


Wstart = 0.01;
Wend = 3.5
ds = 0.025;

i = 2;

Ubws = CONTINUE(@(Xw) SDOF_NLHYST_HBRESFUN(Xw, m, c, k, Fl*FF(i),
					   @(t,x,z,xd) eldryfric(t,x,z,xd,muN0,kt), h, Nt),
		U0, Wstart, Wend, ds, Copt);

figure(1)
plot(Ubws(end,:), sqrt(sum(Ubws(2:3,:).^2, 1))/FF(i), 'o-', 'LineWidth', 2)
