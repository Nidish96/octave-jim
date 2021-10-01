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
dw = 0.01;

Copt = struct('Nmax', 1000, 'dsmin', 0.001);
Dscale = [1e-2*ones(Nhc*GM.Ndofs,1); Wsp(1)];
Copt.Dscale = Dscale;

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
    UCs{fi} = CONTINUE(@(Uw) GM.HBRESFUN(Uw, Fl*fa, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);
%     UCs{fi} = CONTINUE(@(Uw) GM.HBRESFUN(Uw, Fl*fa, h, Nt, 1e-6), U0, Wen, Wst, dw, Copt);
end
%% Save
save(sprintf('./DATA/Jenkins_FRESP_cfg%d.mat',cfg), 'UCs', 'Fl', 'Fas', 'Wst', 'Wen', 'h', 'Nhc');

%% Load
load(sprintf('./DATA/Jenkins_FRESP_cfg%d.mat',cfg), 'UCs', 'Fl', 'Fas', 'Wst', 'Wen');
load(sprintf('./DATA/Jenkins_EPMC_cfg%d.mat',cfg), 'UwxC');

% load(sprintf('./DATA/Jenkins_FRESP.mat',cfg), 'UCs', 'Fl', 'Fas', 'Wst', 'Wen', 'h', 'Nhc');
% load(sprintf('./DATA/Jenkins_EPMC.mat',cfg), 'UwxC');
% %% Plotting
figure(40)
clf()
% plot(UwxC(end-2,:), (10.^UwxC(end,:)).*sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UwxC(1:end-3, :)).^2)), ...
%     '-', 'LineWidth', 2); hold on
aa = gobjects(size(Fas));
for fi=1:length(Fas)
    aa(fi) = plot(UCs{fi}(end,:), ...
        sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UCs{fi}(1:end-1, :)).^2))/Fas(fi), ...
        '.-');
    hold on
    legend(aa(fi), sprintf('F = %.2f MN', Fas(fi)*1e-6));
end
legend(aa(1:end), 'Location', 'northwest')

set(gca, 'yscale', 'linear')
xlim([Wst Wen])
% ylim([0 1e-6])
xlabel('Frequency (rad/s)')
ylabel('RMS Amplitude')