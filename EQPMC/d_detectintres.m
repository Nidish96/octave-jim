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

analyze = true;
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
Nhmax = 1;  % Number of harmonics
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

%% Nonlinear System to Solve
% 1. Parameterize ws as a function of As through EQPMC
% 2. Require n*w1+m*w2 = 0 with n,m being set by user
% 3. (choice 1) given a1, FIND a2 where the above condition is met
% 4. (choice 2) find a1, a2 where the above condition is met. Need
% 	to add an extra condition here to set it up as two equations in
% 	two variables. 
% 	Perhaps we can set MAC between the two mode-shapes to
% 		unity ? Need to try.

%% Forcing
ws = [0.8; pi];
% rng(1);
% hid = randi(length(h)-1,2,1)
% hid = [1; 2];

hid = [find(h(:,1)==0 & h(:,2)==1); find(h(:,1)==1 & h(:,2)==0)]-1;
hfrc = h(1+hid, :);

amps = 1.0*ones(size(hid));  % 20
fext = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

%% Setup QPHB
Nt = 64;

%% Check Gradients w.r.t ais - WRONG; Need to fix!
Fls = zeros(Nhc*2, 2);
Fls(2+(hid(1)-1)*4+1, 1) = 1;
Fls(2+(hid(2)-1)*4+1, 2) = 1;

% rng(1);
% Xv = [X; Wr; xis; -2];
Xv = rand(Nhc*2+2*Nc,1);
Ap = 1.0;
[R, dRdU, dRdlA, dRdash] = MDL.EQPMCRESFUN([Xv; log10(Ap)],  [1; 1], Fls, h, Nt, eps);

hm = 1e-6;
Rp = MDL.EQPMCRESFUN([Xv; log10(Ap)],  [1+hm; 1], Fls, h, Nt, eps);
Rm = MDL.EQPMCRESFUN([Xv; log10(Ap)],  [1-hm; 1], Fls, h, Nt, eps);

(Rp-Rm)/2/hm

Rp = MDL.EQPMCRESFUN([Xv; log10(Ap)],  [1; 1+hm], Fls, h, Nt, eps);
Rm = MDL.EQPMCRESFUN([Xv; log10(Ap)],  [1; 1-hm], Fls, h, Nt, eps);
(Rp-Rm)/2/hm
