clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

%% Setup system (1D Bar)
ell = 1.0;
A = pi*1e-4;
E = 210e9;
rho = 7800;
Ne = 8;

Le = ell/Ne;
Me = 1/3*[2 1;1 2]*Le/2*rho*A;
Ke = 1/2*[1 -1;-1 1]*2/Le*E*A;

M = zeros(Ne+1);
K = zeros(Ne+1);
for e=1:Ne
    M(e:e+1,e:e+1) = M(e:e+1,e:e+1) + Me;
    K(e:e+1,e:e+1) = K(e:e+1,e:e+1) + Ke;
end
Lb = eye(Ne+1);
Lb(:, 1) = [];  % Boundary Condition (DOF 1 fixed)

Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;

%% Damping Factors
Zetas = [0.4; 0.4]*1e-2;

[V, Wsp] = eig(Kb, Mb);
[Wsp, si] = sort(sqrt(diag(Wsp)));
V = V(:, si);
V = V./diag(sqrt(V'*Mb*V))';

if Ne>2
    ab = [1./(2*Wsp(1:length(Zetas))) Wsp(1:length(Zetas))/2]\Zetas;
    Cb = ab(1)*Mb + ab(2)*Kb;
else
    Cb = 2*Zetas(1)/Wsp(1)*Kb;
end

%% Setup model
GM = MDOFGEN(Mb, Kb, Cb, Lb);

% kc = 1e6;
% fnl = @(t,u,ud) deal(kc*u.^3, 3*kc*u.^2, zeros(size(u)));
cfg = 2;

if cfg==1
    % Shows 1:3 mode coupling (cfg=1)
    kt  = 6e6;
    muN = 7.5e5;
elseif cfg==2
    % Doesn't show mode coupling (cfg=2)
    kt  = 2.5e6;
    muN = 7.5e5;
else
    error('unknown cfg')
end

fnl = @(t, u, varargin) JENKFORCE(t, u, kt, muN, varargin{:});
GM = GM.SETNLFUN(2+3, Lb(end,:), fnl);

%% Backbone
h = 0:5;
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^9;

Fl = kron([0; 1; 0; zeros(Nhc-3,1)], Lb(end,:)');

As = -2;
Ae = 2;
da = 0.05;

mds = [1:2];
Copt = struct('Nmax', 1000, 'dsmax', 0.05, 'dsmin', 0.001);
UwxC = cell(size(mds));
dUwxC = cell(size(mds));
for mi=mds
    Uwx0 = [kron([0; 0; 1; zeros(Nhc-3,1)], V(:,mi)); Wsp(mi); 2*Zetas(mi)*Wsp(mi)];
    Dscale = [ones(Nhc*GM.Ndofs,1); Wsp(mi); 2*Zetas(mi)*Wsp(mi); 1.0];
    Copt.Dscale = Dscale;

    [UwxC{mi}, dUwxC{mi}] = CONTINUE(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt, 1e-6), Uwx0, As, Ae, da, Copt);
end
% %% Jaccheck
% opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', false);
% Us = fsolve(@(uwx) GM.EPMCRESFUN([uwx; -0.1], Fl, h, Nt, 1e-6), Uwx0, opt);

%% Save
save(sprintf('./DATA/Jenkins_EPMC_cfg%d.mat',cfg), 'UwxC', 'dUwxC', 'h', 'Nhc', 'mds');

%% Load
load(sprintf('./DATA/Jenkins_EPMC_cfg%d.mat',cfg), 'UwxC', 'mds');

%% Plot
figure(1);
clf();
for mi=mds
    plot(10.^UwxC{mi}(end,:), UwxC{mi}(end-2,:), '.-'); hold on
end

set(gca, 'XScale', 'log')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (rad/s)')

figure(2);
clf()
for mi=mds
    plot(10.^UwxC{mi}(end,:), UwxC{mi}(end-1,:)./(2*UwxC{mi}(end-2,:)), '.-'); hold on
end
set(gca, 'XScale', 'log')
xlabel('Modal Amplitude')
ylabel('Damping Factor')
%%
figure(3)
clf();
aa = gobjects(h(end), 1);
mi = 1;
for hi=1:h(end)
    hv = zeros(Nhc,1);
    hv((h(1)==0)+(hi-1)*2+(1:2)) = 1;
        
	aa(hi) = semilogy(UwxC{mi}(end-2,:), ...
        sqrt(sum((kron(diag(hv),Lb(end,:))*UwxC{1}(1:end-3,:)).^2)), ...
        'LineWidth', 2); hold on
    legend(aa(hi), sprintf('Harmonic %d', hi))
%     plot(UwxC{mi}(end-2,:), sum((kron(diag(hv),Lb(end,:))*UwxC{1}(1:end-3,:)).^2))
end
legend(aa(1:end), 'Location', 'best')
ylim([1e-4 1e0])
%% 
figure(4)
clf()
aa = gobjects(size(mds));
for mi=mds
    aa(mi) = plot((0:Ne)*Le, Lb*V(:,mi), 'o-'); hold on
    legend(aa(mi), sprintf('Mode %d', mi))
end
legend(aa(1:end), 'Location', 'best')
grid on
xlabel('X Coordinate (m)')
ylabel('Deflection (m)')
%% Modal interaction?
figure(10)
clf()
semilogx(10.^UwxC{1}(end,:), UwxC{1}(end-2,:), '.-'); hold on
semilogx(10.^UwxC{2}(end,:), UwxC{2}(end-2,:)/3, '.-')

legend('Mode 1', 'Mode 2 (freq/3)')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (rad/s)')