clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')
addpath('../ROUTINES/export_fig/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

analyze = false;
plotout = false;
%% Parameters
M = eye(2);
K = [2 -1;-1 2];
C = 0.01*K;
b = 0.5;

MDL = MDOFGEN(M, K, C, eye(2));
MDL = MDL.SETNLFUN(1+3, [1 0], @(t, u, ud) deal(b*u.^3, 3*b*u.^2, ud*0));

[V, Wr] = eig(K, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

Nc = 2;  % Number of components
Nhmax = 5;  % Number of harmonics
%% Harmonic Selection
hall = cell(1, Nc);
[hall{:}] = ndgrid(-Nhmax:Nhmax);
hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
h = hall(sum(abs(hall),2)<=Nhmax & sum(hall,2)>=0,:);
% h = hall(abs(sum(hall,2))<=Nhmax & sum(hall,2)>=0,:);

h(sum(h,2)==0 & h(:,1)<=0, :) = [];
h = [zeros(1,Nc); h];
Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));

figure(1)
clf()
plot(hall(:,1), hall(:,2), 'ko', 'MarkerFaceColor', 'w'); hold on
plot(h(:,1), h(:,2), '*'); hold on
grid on
axis equal
legend('All Harmonics', 'Selected Harmonic', 'Location', 'northoutside')

%% Forcing
ws = [0.8; pi];
% rng(1);
% hid = randi(length(h)-1,2,1)
% hid = [1; 2];

hid = [find(h(:,1)==0 & h(:,2)==1); find(h(:,1)==1 & h(:,2)==0)]-1;
hfrc = h(1+hid, :);

amps = 1.0*ones(size(hid));  % 20
fext = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

%% QP HB Simulation
Nt = 64;

[E, dEdw] = QPHARMONICSTIFFNESS(MDL.M, MDL.C, MDL.K, ws, h);
Fl = zeros(Nhc*2, 1);
Fl(2+(hid-1)*4+1) = amps;

X0 = E\Fl;

opt = struct('Display', true, 'ITMAX', 100, 'crit', 30);
% X = NSOLVE(@(U) MDL.QPHBRESFUN([U; 1], ws, Fl, h, Nt, eps), X0, opt);

% fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% X = fsolve(@(U) MDL.QPHBRESFUN([U; 1], ws, Fl, h, Nt, eps), X0, fopt);
%% EQPMC Simulations
Fls = zeros(Nhc*2, 2);
Fls(2+(hid(1)-1)*4+1, 1) = 1;
Fls(2+(hid(2)-1)*4+1, 2) = 1;

xis = [M(:) K(:)]\C(:);  % Fit Proportional Constants (2 term Caughy series)

X0 = zeros(Nhc*2, 1);
X0(2+(hid(1)-1)*4+(1:4)) = kron(V(:,1), [0;1]);
X0(2+(hid(2)-1)*4+(1:4)) = kron(V(:,2), [0;1]);
Xv = [X0; Wr; xis; -2];
% [R, dRdU, dRdlA] = MDL.EQPMCRESFUN(Xv,  [1; 1], Fls, h, Nt, eps); 
Uwx0 = NSOLVE(@(X) MDL.EQPMCRESFUN([X; Xv(end)],  [1; 1], Fls, h, Nt, eps), Xv(1:end-1), opt);

% fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% Uwx = fsolve(@(X) MDL.EQPMCRESFUN([X; Xv(end)],  [1; 1], Fls, h, Nt, eps), Xv(1:end-1), fopt);

%% Continuation
Copt = struct('Nmax', 200, 'Display', 1, 'DynDscale', 0, 'solverchoice', 2, 'angopt', 1e-1);
Copt.Dscale = [kron([1e-2; 1e-1*ones(Nhc-1,1)], ones(2,1)); Wr; 1e-1*ones(2,1); 1.0];
% Copt.Dscale = [kron([1e-4; ones(Nhc-1,1)], ones(2,1)); 1e-2*ones(4,1); 1];

Astart = -2;
Aend = 2;
da = 0.05;
Copt.dsmax = 0.4;

if analyze
    thetas = linspace(0, pi/2, 8); thetas = thetas(2:end-1);
    UwxL = cell(size(thetas));
    for ti=1:length(thetas)
        UwxL{ti} = CONTINUE(@(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [cos(thetas(ti)); sin(thetas(ti))], Fls, h, Nt, eps), ...
            Uwx0, Astart, Aend, da, Copt);
    end
    save('DATA/EQPSURF.mat', 'thetas', 'UwxL')
else
    load('DATA/EQPSURF.mat', 'thetas', 'UwxL')
end

%% Plot 2D Surface

figure(1000)
clf()
for i=1:2
    subplot(2,1,i)
    for ti=1:length(thetas)
        plot3((10.^UwxL{ti}(end,:))*cos(thetas(ti)), (10.^UwxL{ti}(end,:))*sin(thetas(ti)), UwxL{ti}(end-5+i,:), '.-', 'LineWidth', 2); hold on
    end
    set(gca, 'XScale', 'log', 'YScale', 'log')
    grid on
    xlabel('Modal Amp $Q_1$')
    ylabel('Modal Amp $Q_2$')
    zlabel(sprintf('Mode %d Freq', i))
end