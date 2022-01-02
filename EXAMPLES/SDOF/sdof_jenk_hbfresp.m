clc
clear all
addpath('../../ROUTINES/')
addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/SOLVERS/')
addpath('../../ROUTINES/CONTACTMODELS/')  % Find "JENKFORCE.m" here

%% Parameters
mul = 10;

m = 1.0*mul;
k = 4.0*mul;
c = 2*0.005*sqrt(k/m)*mul;

kt = 5.0*mul;
muN = 0.1*mul;

fnl = @(t, u, varargin) JENKFORCE(t, u, kt, muN, varargin{:});

%% Setup Model
GM = MDOFGEN(m, k ,c, 1.0);
GM = GM.SETNLFUN(2+3, 1.0, fnl);

%% HBM Setup Fresp
h = 0:5; 
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^7;
Fl = [0; 1; 0; zeros(Nhc-3,1)];

Wst = sqrt(k/m)*0.5;
Wen = sqrt(k/m)*2;
dw = 0.01;

Copt = struct('Nmax', 1000, 'DynScale', 0);  % Dynamic Scaling Improves solver tolerance but sometimes makes it too difficult to converge unnecessarily

%% Calculate and Plot
figure(2)
clf()
for famp = [0.005 0.01 0.05 0.10 0.15]*mul
    U0 = HARMONICSTIFFNESS(GM.M, GM.C, GM.K, Wst, h)\(Fl*famp);

    GM.Rsc = max(abs(U0), min(abs(U0(U0~=0)))); 
    UCs = PRECOCONT(@(Uw) GM.HBRESFUN(Uw, Fl*famp, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);
    
%     %% NLVib Continuation (If above doesn't work use this)
%     Sopt = struct('jac', 'full', 'dynamicDscale', 1);
%     UCs = solve_and_continue(U0, @(Uw) GM.HBRESFUN(Uw, Fl*famp, h, Nt, 1e-6, 1), Wst, Wen, dw, Sopt);

    %% Plot
    Urms = sqrt([1, 0.5*ones(1,Nhc-1)]*UCs(1:end-1,:).^2);  % RMS Amplitude from Parseval's Thm.
    Uamp1 = abs(UCs(2,:)-1j*UCs(3,:));  % First Harmonic Amplitude
    Uph1 = rad2deg(unwrap(angle(UCs(2,:)-1j*UCs(3,:))));  % First Harmonic Phase

    figure(2)
    subplot(2,2,1)
    plot(UCs(end,:), Uamp1/famp, '-'); hold on
    ylabel('Displacement/Force')

    subplot(2,2,2)
    semilogx(Urms, UCs(end,:), '-'); hold on 
    
    subplot(2,2,3)
    plot(UCs(end,:), Uph1, '-'); hold on
    ylabel('Phase (degs)')
    xlabel('Frequency (rad/s)')
end

%% EPMC Solution
As = -1.8;
Ae = 0.6;
da = 0.01;

% Copt.solverchoice = 2; 
% Copt.arclengthparm = 'K0NormalizedArcLength';

Uwx0 = zeros(Nhc, 1);
Uwx0(2:3) = 1/sqrt(m*mul);  % "Mode" Initialization
Uwx0 = [Uwx0; sqrt(k/m); c/mul];
[UwxC, dUwxC] = CONTINUE(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt, 1e-6), Uwx0, As, Ae, da, Copt);

%% Plot EPMC Solution Over Freps
Ubb = (10.^UwxC(end, :)).*UwxC(1:end-3,:);
Udbb = (HARMONICSTIFFNESS(0, 1, 0, 1, h)*Ubb).*UwxC(end-2,:);  % Harmonics of Velocity

Urms = sqrt([1, 0.5*ones(1,Nhc-1)]*Ubb.^2);  % RMS Amplitude from Parseval's Thm.
Uamp1 = abs(Ubb(2,:)-1j*Ubb(3,:));  % First Harmonic Amplitude
Uph1 = rad2deg(angle(Ubb(2,:)-1j*Ubb(3,:)));  % First Harmonic Phase

Feff = (kron(eye(Nhc), GM.M)*Udbb).*UwxC(end-1,:);
Feff1 = abs(Feff(2,:)-1j*Feff(3,:));  % First Harmonic Amplitude
Feffrms = sqrt([1, 0.5*ones(1,Nhc-1)]*Feff.^2);  % RMS Amplitude from Parseval's Thm.

figure(2)
subplot(2,2,1)
hold on;
plot(UwxC(end-2,:), Uamp1./(Feffrms*sqrt(2)), 'k-', 'LineWidth', 1)
subplot(2,2,1); set(gca, 'yscale', 'log')
title('Frequency Response')

subplot(2,2,3)
hold on;
plot(UwxC(end-2,:), Uph1, 'k-', 'LineWidth', 1.2)
title('Frequency Reponse Phase')

subplot(2,2,2)
hold on
semilogx(Urms, UwxC(end-2,:), 'k-', 'LineWidth', 1.2);
ylabel('Natural Frequency')
title('Frequency Backbone')

subplot(2,2,4)
semilogx(Urms, UwxC(end-1,:)./(2*UwxC(end-2,:)), 'k-', 'LineWidth', 1.2);
xlabel('RMS Amplitude')
ylabel('Effective Damping')
title('Damping Backbone')