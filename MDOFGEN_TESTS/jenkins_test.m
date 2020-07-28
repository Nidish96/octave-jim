clc
clear all 

%%
kt = 1.0;
muN = 0.5;

h = 0:1;
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^3;
t = linspace(0, 2*pi, Nt+1)';  t(end) = [];

fun = @(t, u, varargin) JENKFORCE(t, u, kt, muN, varargin{:});

u = cos(t);
f = zeros(size(t));
df = zeros(size(t), 1, 1, Nhc);
f1 = zeros(size(t));
df1 = zeros(size(t), 1, 1, Nhc);
f2 = zeros(size(t));
df2 = zeros(size(t), 1, 1, 1);
f3 = zeros(size(t));
df3 = zeros(size(t), 1, 1, 1);
for ti=1:Nt
  tm1 = mod(ti-2, Nt)+1;
  
  [f(ti), df(ti, :, :, :)] = JENKFORCE(t(ti), u(ti), kt, muN, h, t(tm1), u(tm1), f(tm1), df(tm1, :, :, :));
  [f1(ti), df1(ti, :, :, :)] = fun(t(ti), u(ti), h, t(tm1), u(tm1), f(tm1), df(tm1, :, :, :));
  
  [f2(ti), df2(ti, :, :, :)] = JENKFORCE(t(ti), u(ti), kt, muN);
  [f3(ti), df3(ti, :, :, :)] = fun(t(ti), u(ti));
end

figure(1)
clf()
plot(t, f, 'k-'); hold on 
plot(t, f1, 'o');
plot(t, f2, 'b-');
plot(t, f3, 'o');
