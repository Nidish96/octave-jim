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

analyze = true;
plotout = false;
%% Parameters
M = eye(2);
K = [2 -1;-1 2];
C = 0.01*K;

MDL = MDOFGEN(M, K, C, eye(2));

kt = 4;
muN = 2;
MDL = MDL.SETNLFUN(2+3, [1 0], @(t, u, varargin) JENKNL(t, u, kt, muN, varargin{:}), [], 4);
K0 = K+kt*[1 0;0 0]; % Linearized stiffness
[V, Wr] = eig(K0, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

Nc = 2;  % Number of components
Nhmax = 1;  % Number of harmonics
%% Harmonic Selection
h = HSEL(Nhmax, [1 2]);
Nhc = sum(all(h==0,2) + 2*any(h~=0,2));

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

[E, dEdw] = QPHARMONICSTIFFNESS(MDL.M, MDL.C, K0, ws, h);
Fl = zeros(Nhc*2, 1);
Fl(2+(hid-1)*4+1) = amps;

X0 = E\Fl;

% opt = struct('Display', true, 'ITMAX', 100, 'crit', 30);
% X = NSOLVE(@(U) MDL.QPHBRESFUN([U; ws(:)], Fl, h, Nt, eps), X0, opt);

% fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% X = fsolve(@(U) MDL.QPHBRESFUN([U; ws(:)], Fl, h, Nt, eps), X0*0, fopt);

%% EQPMC Simulations
xis = [M(:) K0(:)]\C(:);  % Fit Proportional Constants (2 term Caughy series)

[rinds0,zinds,hinds,rinds,iinds] = HINDS(2, h);
Fls = zeros(Nhc*2, 2);
Fls(rinds(1), 1) = 1;
Fls(rinds(3), 2) = 1;

X0 = zeros(Nhc*2, 1);
X0([rinds(1:2) iinds(1:2)]) = kron([0;1], V(:,1));
X0([rinds(3:4) iinds(3:4)]) = kron([0;1], V(:,2));
Xv = [X0; Wr; xis; -1];
% [R, dRdU, dRdlA] = MDL.EQPMCRESFUN(Xv,  [1; 1], Fls, h, Nt, eps); 

% opt = struct('Display', true, 'ITMAX', 100, 'lsrch', 1);
% Uwx = NSOLVE(@(X) MDL.EQPMCRESFUN([X; Xv(end)],  [1; 1], Fls, h, Nt, eps), Xv(1:end-1), opt);

fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display', 'iter');
Uwx = fsolve(@(X) MDL.EQPMCRESFUN([X; Xv(end)],  [1; 1], Fls, h, Nt, eps), Xv(1:end-1), fopt);

%% Continuation
Copt = struct('Nmax', 1000, 'Display', 1);
Astart = -2;
Aend = 2;
da = 0.1;

Uwx0 = Xv(1:end-1);
UwxL = CONTINUE(@(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [1; 1], Fls, h, Nt, eps), Uwx0, Astart, Aend, da, Copt);

%% 
figure(2)
clf()
loglog(10.^UwxL(end,:), UwxL(end-4:end-3,:), '.-')
xlabel('Amplitude Scaling')
ylabel('Frequency (rad/s)')

figure(3)
clf()
mamps = [diag(UwxL(2+(hid(1)-1)*4+(1:2),:)'*MDL.M*UwxL(2+(hid(1)-1)*4+(1:2),:))' + diag(UwxL(2+(hid(1)-1)*4+(3:4),:)'*MDL.M*UwxL(2+(hid(1)-1)*4+(3:4),:))';
    diag(UwxL(2+(hid(2)-1)*4+(1:2),:)'*MDL.M*UwxL(2+(hid(2)-1)*4+(1:2),:))' + diag(UwxL(2+(hid(2)-1)*4+(3:4),:)'*MDL.M*UwxL(2+(hid(2)-1)*4+(3:4),:))'];
kamps = [diag(UwxL(2+(hid(1)-1)*4+(1:2),:)'*MDL.K*UwxL(2+(hid(1)-1)*4+(1:2),:))' + diag(UwxL(2+(hid(1)-1)*4+(3:4),:)'*MDL.K*UwxL(2+(hid(1)-1)*4+(3:4),:))';
    diag(UwxL(2+(hid(2)-1)*4+(1:2),:)'*MDL.K*UwxL(2+(hid(2)-1)*4+(1:2),:))' + diag(UwxL(2+(hid(2)-1)*4+(3:4),:)'*MDL.K*UwxL(2+(hid(2)-1)*4+(3:4),:))'];
semilogx(10.^UwxL(end,:), UwxL(end-2,:).*mamps+UwxL(end-1,:).*kamps, '.-')
xlabel('Amplitude Scaling')
ylabel('Damping Coefficients')

% %% Check Gradients
% rng(1);
% Xv = [X; Wr; xis; -2];
% Xv = rand(Nhc*2+2*Nc+1,1);
% [R, dRdU, dRdlA] = MDL.EQPMCRESFUN(Xv,  [1; 1], Fls, h, Nt, eps); 
% 
% hm = 1e-9;
% hv = zeros(Nhc*MDL.Ndofs+2*Nc+1, 1);
% Jnum = zeros(Nhc*MDL.Ndofs+2*Nc, Nhc*MDL.Ndofs+2*Nc+1);
% for hi=1:(Nhc*MDL.Ndofs+2*Nc+1)
%     hv(hi) = hm;
%     Rp = MDL.EQPMCRESFUN(Xv+hv,  [1; 1], Fls, h, Nt, eps); 
%     Rm = MDL.EQPMCRESFUN(Xv-hv,  [1; 1], Fls, h, Nt, eps); 
%     hv(hi) = 0;
%     
%     Jnum(:, hi) = (Rp-Rm)/2/hm;
% end
% disp([max(max(abs((Jnum(:,1:end-1)-dRdU)))) max(abs(Jnum(:,end)-dRdlA))])
