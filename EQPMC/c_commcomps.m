clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')
addpath('../ROUTINES/export_fig/')
% set(0,'defaultAxesTickLabelInterpreter', 'default');
% set(0,'defaultTextInterpreter','latex');
% set(0, 'DefaultLegendInterpreter', 'latex');
% set(0,'defaultAxesFontSize',13);

%% Parameters
M = eye(2);
K = [2 -1;-1 2];
C = 0.01*K;
b = 0.5;

MDL = MDOFGEN(M, K, C, eye(2));
% MDL = MDL.SETNLFUN(1+3, [1 0], @(t, u, ud) deal(b*u.^3, 3*b*u.^2, ud*0));

[V, Wr] = eig(K, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

Nc = 2;  % Number of components
Nhmax = 3;  % Number of harmonics
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
ws = [0.25; 0.25];
% rng(1);
% hid = randi(length(h)-1,2,1)
% hid = [1; 2];

hid = [find(h(:,1)==1 & h(:,2)==0); find(h(:,1)==0 & h(:,2)==1)]-1;
hfrc = h(1+hid, :);

amps = 1.0*ones(size(hid));  % 20
fext = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

%% QP HB Simulation
Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
Nt = 64;

E = QPHARMONICSTIFFNESS(MDL.M, MDL.C, MDL.K, ws, h);
Fl = zeros(Nhc*2, 1);
Fl(2+(hid-1)*4+1) = amps;

X0 = E\Fl;

opt = struct('Display', true, 'ITMAX', 100, 'crit', 30);
X = NSOLVE(@(U) MDL.QPHBRESFUN([U; 1], ws, Fl, h, Nt, eps), X0, opt);

%% Regular HB Simulation
hh = 0:3; Nhch = sum((hh==1)+2*(hh~=0));

Eh = HARMONICSTIFFNESS(MDL.M, MDL.C, MDL.K, ws(1), hh);
Flh = zeros(Nhch*2, 1);
Flh(2+(1-1)*4+1) = Flh(2+(1-1)*4+1)+amps(1);
Flh(2+(1-1)*4+1) = Flh(2+(1-1)*4+1)+amps(2);

X0h = Eh\Flh;
Xh = NSOLVE(@(U) MDL.HBRESFUN([U; ws(1)], Flh, hh, Nt, eps), X0h, opt);
