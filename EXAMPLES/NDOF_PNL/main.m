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
% Nyquist is around 420 Hz

C = 0.68*M + 3.75e-6*K;
% diag(V'*C*V)./(2*sqrt(diag(V'*K*V)))

%% Setup Forcing
Lf = zeros(Ndofs, 1);
Lf(end) = 1;
Fex = @(t) Lf.*cos(5*2*pi*t)*20;

%% Setup Nonlinearity 
p = 3;  % Highest order of nonlinearity
qq = dec2base((0:(p+1)^Ndofs-1)', p+1);
qq = arrayfun(@(c) str2num(c')', qq);

qq = qq(~(sum(qq, 2)<=1 | sum(qq, 2)>p), :);

E = rand(Ndofs, size(qq,1))*10;

Fnl = @(t, u, ud, tp) deal(E*prod(u'.^qq, 2), ...
    E*((prod(u'.^qq, 2)./(u'.^qq+eps)).*(qq.*(u'.^max(qq-1,0)))), ...  % adding eps to denominator of first term to get rid of div by 0 issues
    zeros(Ndofs));

% %% Check Gradients
% rng(2)
% U0 = rand(Ndofs, 1);
% [F0, J0, ~] = Fnl(0, U0, U0*0, 0);
% 
% hi = 1;
% hm = 1e-6;
% hv = zeros(Ndofs, 1);
% hv(hi) = 1;
% [Fp, ~, ~] = Fnl(0, U0+hv*hm, U0*0, 0);
% [Fm, ~, ~] = Fnl(0, U0-hv*hm, U0*0, 0);
% hv(hi) = 0;
% 
% [(Fp-Fm)/(2*hm) J0(:, hi)]

%% 
MDL_lin = MDOFGEN(M, K, C, eye(Ndofs));  % Linear System (no non-linearities

MDL = MDL_lin.SETNLFUN(1+3, eye(Ndofs), Fnl);  % Nonlinear System (along with polynomial terms)

%% Static Equilibrium
U0 = -ones(Ndofs, 1);  % We want a non-trivial equilibrium

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
Ustat = fsolve(@(U) MDL.STATRESFUN(U, U0*0), U0, opt);

%% 
Tmax = 20.0;
fsamp = 1000; % Nyquist is around 420 Hz. We're doing more than twice that. 4-5 times is usually necessary for converged solution

Hopt = struct('alpha', 0, 'beta', 1/4, 'gamma', 1/2, ...
      'Display', 'waitbar');  % Get rid of Display (set to 'off') to improve speed
[~, U, Ud, Udd] = MDL_lin.HHTAMARCH(0, Tmax, 1/fsamp, zeros(Ndofs,1), zeros(Ndofs,1), ...
    Fex, Hopt);
[T, Unl, Udnl, Uddnl] = MDL.HHTAMARCH(0, Tmax, 1/fsamp, zeros(Ndofs,1), zeros(Ndofs,1), ...
    Fex, Hopt);

%%
figure(1)
clf()
subplot(2,1,1)
plot(T, Lf'*U, '-', 'LineWidth', 2); hold on
plot(T, Lf'*Unl, '-', 'LineWidth', 2)
plot(T, ones(size(T))*(Lf'*Ustat), 'k-', 'LineWidth', 2)
ylabel('Displacement (m)')
subplot(2,2,3)
plot(T, Lf'*Fex(T), '.-')
xlabel('Time (s)')
ylabel('Force (N)')
subplot(2,2,4)
plot(Lf'*U, Lf'*Ud, Lf'*Unl, Lf'*Udnl); hold on
plot([Lf'*Ustat; 0], [0; 0], 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 2)
legend('Linear System', 'Nonlinear System', 'Nonlinear Equilibria')
xlabel('Displacement (m)')
ylabel('Velocity (m/s)')