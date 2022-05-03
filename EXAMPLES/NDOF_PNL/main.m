clc
clear all
addpath('../../ROUTINES/')

%% 
Ndofs = 6;

rng(1);
M = diag(rand(Ndofs,1));

K = zeros(Ndofs);
for e=1:Ndofs-1
    K(e:e+1, e:e+1) = K(e:e+1, e:e+1) + [1 -1;-1 1]*100;
end
K(1,1) = 2*K(1,1);

C = 0.68*M + 3.75e-6*K;
% diag(V'*C*V)./(2*sqrt(diag(V'*K*V)))

%% Setup Forcing
Lf = zeros(Ndofs, 1);
Lf(end) = 1;
Fex = @(t) Lf.*cos(5*2*pi*t)*100;

%% Setup Nonlinearity 
p = 3;  % Highest order of nonlinearity
qq = dec2base((0:4^Ndofs-1)', p+1);
qq = arrayfun(@(c) str2num(c')', qq);

qq = qq(~(sum(qq, 2)<=1 | sum(qq, 2)>3), :);

E = rand(Ndofs, size(qq,1));

Fnl = @(t, u, ud, tp) deal(E*prod(u'.^qq, 2), ...
    E*(qq.*(u'.^max(qq-1,0))), ...
    zeros(Ndofs));

%% 
MDL = MDOFGEN(M, K, C, eye(Ndofs));

MDL = MDL.SETNLFUN(1+3, eye(Ndofs), Fnl);

%% 
Tmax = 20.0;
fsamp = 100;

Hopt = struct('alpha', 0, 'beta', 1/4, 'gamma', 1/2, ...
      'Display', 'waitbar');
[T, U, Ud, Udd] = MDL.HHTAMARCH(0, Tmax, 1/fsamp, zeros(Ndofs,1), zeros(Ndofs,1), ...
    Fex, Hopt);

%%
figure(1)
clf()
plot(T, Lf'*U, '.-')
yyaxis right
plot(T, Fex(T), '-')