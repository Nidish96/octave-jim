clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC')
addpath('../ROUTINES/SOLVERS')

%% Description
% Using a 2dof model to determine if the implemented functions for ELDRYFRICT2D provide correct gradients
% dof 1 -> tangential; dof 2 -> normal

%% Model
M = [1 0;0 1];
K = [4 0;0 4];
C = [2*0.005*2 0;0 2*0.005*2];
    
gap = 0.0;
mu = 0.9;
kt = 5;
kn = 50;

MDL = MDOFGEN(M, K, C, eye(2));

fnl = @(t, u, varargin) ELDRYFRICT2D(t, u, kt, kn, mu, gap, varargin{:});
MDL = MDL.SETNLFUN(2+5, eye(2), fnl, eye(2));

%% Static Analysis
Fn = [0; 1];
Ft = [1; 0];

% Static forces
fns = 10;
fts = 0;

U0 = MDL.K\(Fn*fns+Ft*fts);

opts = struct('reletol', 1e-6, 'Display', true);
[Ustat, ~, ~, ~, J0] = NSOLVE(@(U) MDL.STATRESFUN(U, Fn*fns+Ft*fts), U0, opts);

%% HBM
h = [0 1 2 3];
Nhc = sum((h==0)+2*(h~=0));

Nd = 2;
Nt = 2^9;

Wst = 0.1;
Wen = 5;

% Wst = 5;
% Wen = 0.1;

dw = 0.1;

% Dynamic forces
fnd = 0.5;  % 2, 0.5
ftd = 0.5;  % 0.1 ,0.5, 1.0, 2.0, 4

Fl = [Fn*fns+Ft*fts;
    Fn*fnd+Ft*ftd;
    zeros(Nd*(Nhc-2),1)];

U0 = HARMONICSTIFFNESS(MDL.M,MDL.C,J0,Wst,h)\Fl;

Copt = struct('Nmax', 5000, 'itDisplay', false, 'angopt', 1e-1,...
    'arclengthparm', 'orthogonal', 'lsrch', 2);
% Copt.Dscale = [ones(Nd*Nhc,1)*5e0; 1];
% UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

UC = CONTINUE(@(Uw) TMPHBRESFUN(Uw, M, C, K, Fl, kt, kn, mu, Nt, h), U0, Wst, Wen, dw, Copt);

% Sopt = struct('stepmax',1e3, 'jac', 'off');
% UC = solve_and_continue(U0, @(Uw) TMPHBRESFUN(Uw, M, C, K, Fl, kt, kn, mu, Nt, h), Wst, Wen, dw, Sopt);

% Nw = 400;
% UC = zeros(Nd*Nhc+1, Nw);
% UC(end,:) = linspace(Wst, Wen, Nw);
% opts = optimset('Jacobian', 'on', 'Display', 'off');
% u0 = HARMONICSTIFFNESS(MDL.M,MDL.C,J0,UC(end,1),h)\Fl;
% 
% for iw=1:Nw
% %     [UC(1:end-1,iw), ~, eflg, ~, Jac] = fsolve(@(u) MDL.HBRESFUN([u; UC(end,iw)], Fl, h, Nt, 1e-6), u0, opts);
%     [UC(1:end-1,iw), ~, eflg, ~, Jac] = fsolve(@(u) TMPHBRESFUN([u; UC(end,iw)], M, C, K, Fl, kt, kn, mu, Nt, h), u0, opts);
%     
%     if eflg<=0
%         keyboard
%     end
%     
%     [~,Jac,Jw] = MDL.HBRESFUN(UC(:,iw), Fl, h, Nt, 1e-6);
%     if iw<Nw
%         u0 = UC(1:end-1,iw) - (Jac\Jw)*(UC(end,iw+1)-UC(end,iw));
%     end
%     fprintf('%d/%d: %d\n', iw,Nw, eflg);
% end

figure(100)
% clf();
subplot(2,1, 1)
% plot(UC(end,:), 20*log10(abs(UC(3:4,:)-UC(5:6,:)*1j)), '.-')
plot(UC(end,:), (abs(UC(3:4,:)-UC(5:6,:)*1j))/fnd, '.-'); hold on
xlabel('Frequency (Hz)')
ylabel('Response (m)')

subplot(2,1, 2)
plot(UC(end,:), rad2deg(angle(UC(3:4,:)-UC(5:6,:)*1j)), '.-'); hold on
xlabel('Frequency (Hz)')
ylabel('Phase (degs)')

%% Time integration
bw = 20/2/pi;
fampt = 0.01;
fampn = 100;
fex = @(t) sin(2*pi*bw*t).^2.*(t<=1.0/2/bw);

fsamp = 64;  dt = 1/fsamp;
T0 = 0;
T1 = 100;
T = (T0:dt:T1)';

fext = fex(T);

FEXv = Fn.*(fext'*fampn+fns) + Ft.*(fext'*fampt+fts);

opts = struct('Display', 'waitbar');
tic
[T, U, Ud, Udd] = MDL.HHTAMARCH(T0, T1, dt, Ustat, [0;0], ...
    FEXv, opts);
toc

figure(1)
clf()
plot(T, Udd(1,:), '.-'); hold on
plot(T, Udd(2,:)); hold on

legend('X', 'N')