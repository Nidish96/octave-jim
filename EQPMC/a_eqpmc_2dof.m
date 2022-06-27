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
b = 0.5;

MDL = MDOFGEN(M, K, C, eye(2));
MDL = MDL.SETNLFUN(1+3, [1 0], @(t, u, ud) deal(b*u.^3, 3*b*u.^2, ud*0));

[V, Wr] = eig(K, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

Nc = 2;  % Number of components
Nhmax = 7;  % Number of harmonics
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
Copt = struct('Nmax', 1000, 'Display', 1);
Astart = -2;
Aend = 2;
da = 0.1;

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
plot(10.^UwxL(end,:), UwxL(end-2,:).*mamps+UwxL(end-1,:).*kamps, '.-')
xlabel('Amplitude Scaling')
ylabel('Frequency (rad/s)')

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