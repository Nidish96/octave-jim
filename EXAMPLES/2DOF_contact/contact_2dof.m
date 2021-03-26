clc
clear all
addpath('../../ROUTINES/')
addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/SOLVERS/')
addpath('../../ROUTINES/CONTACTMODELS/')  % Find "ELDRYFRICT2D.m" here

%% Parameters
m = 1.0;
kx = 4.0;  % resonance at 2 rad/s
cx = 2*0.05*sqrt(kx/m);
kn = 100.0;  % resonance at 10 rad/s
cn = 2*0.1*sqrt(kn/m);

kt = 5.0;  % stuck resonance at 3 rad/s
kn = 21.0;  % contact resonance at 11 rad/s
mu = 0.85;
gap = 0.0;

M = [m 0;0 m];
K = [kx 0;0 kn];
C = [cx 0;0 cn];

Fstat = [0; 100];  % Static force vector (pressing down in the normal spring
Fvec = [1; 0.1];  % Dynamic force

Ustat = (K+diag([kt; kn]))\Fstat;

%% Create Class Object
GM = MDOFGEN(M, K, C, eye(2));

fnl = @(t, u, varargin) ELDRYFRICT2D(t, u, kt, kn, mu, gap, varargin{:});
GM = GM.SETNLFUN(2+3, eye(2), fnl);

%% HBM
h = 0:3;
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^9;
Wst = sqrt(kx/m)*0.5;
Wen = sqrt(kx/m)*2;
dw = 0.1;

Copt = struct('Nmax', 100, 'DynScale', 1, 'arclengthparm', 'arclength');

Famps = [1 10 50 100 200];
for fi=1:length(Famps)
    Fl = [Fstat; Fvec*Famps(fi); zeros(GM.Ndofs*(Nhc-2),1)];
    U0 = HARMONICSTIFFNESS(GM.M, GM.C, GM.K+diag([kt,kn]), Wst, h)\Fl;  % "stuck, contact" initial guess

    Copt.Dscale = [max(abs(U0), min(abs(U0(U0~=0)))); Wst];

    GM.Rsc = max(abs(U0), min(abs(U0(U0~=0)))); 
    % Use this callback to view plot as it evolves
%     figure(5)
%     Copt.CallbackFun = @(Ulc, dUlds, Ulpred, Copt) plot(Ulc(end,:), abs(Ulc(3,:)+1j*Ulc(5,:)), '.-', Ulc(end, end-1), abs(Ulc(3,end-1)+1j*Ulc(5,end-1)), 'o', Ulpred(end), abs(Ulpred(3)+1j*Ulpred(5)), '*');
%     clf();
    UCs{fi} = PRECOCONT(@(Uw) GM.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

    % Nlvib Continuation (try if PRECOCONT fails) (from Malte Krack & Johann Gross)
    % Sopt = struct('jac', 'full', 'dynamicDscale', 1);
    % UCs = solve_and_continue(U0, @(Uw) GM.HBRESFUN(Uw, Fl, h, Nt, 1e-6, 1), Wst, Wen, dw, Sopt);
end

%% EPMC
As = -1;
Ae = 3;
da = 0.01;

Uwx0 = [Ustat; 0; 0; 1; 0; zeros((Nhc-3)*2,1); sqrt((kx+kt)/m); cx];

Copt.Dscale = [max(abs(Uwx0), min(abs(Uwx0(Uwx0~=0)))); 1];
% Use this callback to view plot as it evolves
% figure(5)
% Copt.CallbackFun = @(Uwxa, dUwxas, Uwxapred, Copt) plot(Uwxa(end,:), Uwxa(end-2,:), '.-', Uwxa(end, end-1), Uwxa(end-2, end-1), 'o', Uwxapred(end), Uwxapred(end-2), '*');  % Live Plot Frequency
% Copt.CallbackFun = @(Uwxa, dUwxas, Uwxapred, Copt) plot(Uwxa(end,:), Uwxa(end-1,:)./(2*Uwxa(end-2,:)), '.-', Uwxa(end, end-1), Uwxa(end-1, end-1)./(2*Uwxa(end-2, end-1)), 'o', Uwxapred(end), Uwxapred(end-1)/(2*Uwxapred(end-2)), '*');  % Live Plot Damping
% clf();
[UwxC, dUwxC] = PRECOCONT(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt, 1e-6), Uwx0, As, Ae, da, Copt);

% Nlvib Continuation (try if PRECOCONT fails) (from Malte Krack & Johann Gross)
% Sopt = struct('jac', 'full', 'dynamicDscale', 1);
% UwxCs = solve_and_continue(Uwx0, @(Uwxa) GM.EPMCRESFUN(Uwxa, Fl, h, Nt, 1e-6, 1), As, Ae, da, Sopt);

%% Post process EPMC results
Asc = kron([ones(1, size(UwxC,2)); ones(Nhc-1,1)*(10.^UwxC(end,:))], ones(GM.Ndofs,1));  % Modal Scaling applied only to harmonics (not to zero harmonic)

Ubb = UwxC(1:end-3, :).*Asc;  % Displacements
Feff = (HARMONICSTIFFNESS(zeros(GM.Ndofs), GM.M, zeros(GM.Ndofs), 1, h)*Ubb).*(UwxC(end-1,:).*UwxC(end-2,:));  % "Effective Modal Forces"

%%
figure(4)
clf()
for fi=1:length(Famps)
    Ux = UCs{fi}(1:2:end-1,:);
    Un = UCs{fi}(2:2:end-1,:);

    figure(4)
    subplot(2,2,1)
    plot(UCs{fi}(end,:), abs(Ux(2,:)+1j*Ux(3,:))/Famps(fi), '-'); hold on
    xlim([Wst Wen])
    
    subplot(2,2,3)
    plot(UCs{fi}(end,:), rad2deg(angle(Ux(2,:)-1j*Ux(3,:))), '-'); hold on
    xlim([Wst Wen])

    subplot(2,2,2)
    semilogx(sqrt([1, 0.5*ones(1,Nhc-1)]*Ux.^2), UCs{fi}(end,:), '-'); hold on
end

Ux = Ubb(1:2:end,:);  Feffx = Feff(1:2:end,:);
Un = Ubb(2:2:end,:);  Feffy = Feff(1:2:end,:);
subplot(2,2, 1)
plot(UwxC(end-2,:), abs((Ux(2,:)-1j*Ux(3,:))./(Feffx(2,:)-1j*Feffx(3,:))), 'k-', 'LineWidth', 1.2)
ylabel('H1 FRF Amplitude')

subplot(2,2, 3)
plot(UwxC(end-2,:), rad2deg(angle((Ux(2,:)-1j*Ux(3,:))./(Feffx(2,:)-1j*Feffx(3,:)))), 'k-', 'LineWidth', 1.2)
ylabel('H1 FRF Phase (degs)')
xlabel('Frequency (rad/s)')

subplot(2,2, 2)
semilogx(sqrt([1, 0.5*ones(1,Nhc-1)]*Ux.^2), UwxC(end-2,:), 'k-', 'LineWidth', 1.2); hold on
xlim(10.^[As Ae])
ylim([Wst Wen])
ylabel('Frequency (rad/s)')

subplot(2,2, 4)
semilogx(sqrt([1, 0.5*ones(1,Nhc-1)]*Ux.^2), UwxC(end-1,:)./(2*UwxC(end-2,:)), 'k-', 'LineWidth', 1.2); hold on
xlim(10.^[As Ae])
ylabel('Effective Damping Factor')
xlabel('RMS Amplitude')