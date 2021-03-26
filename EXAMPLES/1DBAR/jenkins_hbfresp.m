clc
clear all
addpath('../../ROUTINES/')
addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/SOLVERS/')
addpath('../../ROUTINES/CONTACTMODELS/')  % Find "JENKFORCE.m" here

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
cfg = 1;

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

%% HBM Fresp
h = 0:5;
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^9;
Fl = kron([0; 1; 0; zeros(Nhc-3,1)], Lb(end,:)');

Wst = 7500;
Wen = 9000;
dw = 200;

Copt = struct('Nmax', 1000, 'dsmin', 0.001, 'DynScale', 1);
% Dscale = [1e-2*ones(Nhc*GM.Ndofs,1); Wsp(1)];
% Copt.Dscale = Dscale;

% Dscale = [1e0*ones(Nhc*GM.Ndofs,1); Wsp(1)];
% Copt.Dscale = Dscale;
% Copt.dsmax = 1;

Fas = [1e5 5e5 10e5 20e5 35e5];
% UCs = cell(size(Fas));

fa = 25e5;  % 1, 5, 25
for fi=1:length(Fas)
    fa = Fas(fi);
    U0 = HARMONICSTIFFNESS(GM.M, GM.C, GM.K, Wst, h)\(Fl*fa);
    
    if fi==1
        Copt.dsmax = dw;
    else
        Copt.dsmax = dw*5;
    end
%     UCs{fi} = CONTINUE(@(Uw) GM.HBRESFUN(Uw, Fl*fa, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);
%     UCs{fi} = CONTINUE(@(Uw) GM.HBRESFUN(Uw, Fl*fa, h, Nt, 1e-6), U0, Wen, Wst, dw, Copt);
    UCs{fi} = PRECOCONT(@(Uw) GM.HBRESFUN(Uw, Fl*fa, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);
end
%% Save
% save(sprintf('./DATA/Jenkins_FRESP_cfg%d.mat',cfg), 'UCs', 'Fl', 'Fas', 'Wst', 'Wen', 'h', 'Nhc');

%% Load
% load(sprintf('./DATA/Jenkins_FRESP_cfg%d.mat',cfg), 'UCs', 'Fl', 'Fas', 'Wst', 'Wen');
% load(sprintf('./DATA/Jenkins_EPMC_cfg%d.mat',cfg), 'UwxC');

% load(sprintf('./DATA/Jenkins_FRESP.mat',cfg), 'UCs', 'Fl', 'Fas', 'Wst', 'Wen', 'h', 'Nhc');
% load(sprintf('./DATA/Jenkins_EPMC.mat',cfg), 'UwxC');

%% EPMC
As = -2;
Ae = 2;
da = 0.02;

Copt = struct('Nmax', 1000, 'DynScale', 1, 'crit', 6);

Uwx0 = [kron([0; 0; 1; zeros(Nhc-3,1)], V(:,1)); Wsp(1); 2*Zetas(1)*Wsp(1)];
[UwxC, dUwxC] = PRECOCONT(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt, 1e-6), Uwx0, As, Ae, da, Copt);

Ubb = UwxC(1:end-3, :).*(10.^UwxC(end,:));  % Displacements
Feff = (HARMONICSTIFFNESS(zeros(GM.Ndofs), GM.M, zeros(GM.Ndofs), 1, h)*Ubb).*(UwxC(end-1,:).*UwxC(end-2,:));  % "Effective Modal Forces"

%% Plotting
figure(30)
clf()
% plot(UwxC(end-2,:), (10.^UwxC(end,:)).*sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UwxC(1:end-3, :)).^2)), '-', 'LineWidth', 2)
aa = gobjects(length(Fas)+1,1);
for fi=1:length(Fas)
    Urms = sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UCs{fi}(1:end-1, :)).^2));
    Uh1 = kron([0, 1, -1j, zeros(1,Nhc-3)], GM.NLTs.L)*UCs{fi}(1:end-1,:);
    
    subplot(2,2,1)
    aa(fi) = plot(UCs{fi}(end,:), abs(Uh1), '-'); hold on
    legend(aa(fi), sprintf('F = %.2f MN', Fas(fi)*1e-6));
    
    subplot(2,2,2)
    plot(Urms, UCs{fi}(end,:), '-'); hold on
    
    subplot(2,2,3)
    plot(UCs{fi}(end,:), rad2deg(angle(Uh1)), '-'); hold on
end
subplot(2,2,1)
Ua1 = kron([0, 1, 1j, zeros(1, Nhc-3)], GM.NLTs.L)*Ubb;
aa(end)=plot(UwxC(end-2,:), abs(Ua1), 'k-', 'LineWidth', 1.2);
legend(aa(end), 'EPMC');

ll=legend(aa(1:end), 'Location', 'best');
set(ll, 'Visible', 'off')

set(gca, 'yscale', 'log')
xlim([Wst Wen])
ylabel('H1 Displacement')

subplot(2,2,2)
Urms = sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*Ubb).^2));
plot(Urms, UwxC(end-2,:), 'k-');
ylim([Wst/1.05 Wen*1.05])
set(gca, 'xscale', 'log')
ylabel('Natural Frequency')
xlim(10.^[As Ae])

subplot(2,2,3)
plot(UwxC(end-2,:), rad2deg(angle(Ua1)), 'k-', 'LineWidth', 1.2);
ylabel('H1 Phase (degs)')
xlim([Wst Wen])
xlabel('Frequency (rad/s)')

subplot(2,2,4)
plot(Urms, UwxC(end-1,:)./(2*UwxC(end-2,:)), 'k-');
set(gca, 'xscale', 'log')
ylabel('Natural Frequency')
xlabel('RMS Response Amplitude')
xlim(10.^[As Ae])