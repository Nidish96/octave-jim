clc
clear all
addpath('../../ROUTINES/')
addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/SOLVERS/')

%% Parameters
mul = 10;

m = 1.0*mul;
k = 4.0*mul;
c = 2*0.005*sqrt(k/m)*mul;
bt = 0.1*mul;
fnl = @(t,u,ud) deal(bt*u.^3, 3*bt*u.^2, zeros(size(u)));

%% Setup Model
GM = MDOFGEN(m, k ,c, 1.0);
GM = GM.SETNLFUN(1+3, 1.0, fnl);

%% HBM Steady State Forced Response
h = 0:5; 
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^9;
Fl = [0; 1; 0; zeros(Nhc-3,1)];

Wst = sqrt(k/m)*0.5;
Wen = sqrt(k/m)*1.5;
dw = 0.01;

Copt = struct('Nmax', 1000, 'DynScale', 1);  % Dynamic Scaling Improves solver tolerance but sometimes makes it too difficult to converge unnecessarily
Copt.solverchoice = 2; % Test with 1 and 2

%% Calculate and Plot
figure(1)
clf()
for famp = [0.01 0.05 0.10 0.15 0.20]*mul
    U0 = HARMONICSTIFFNESS(GM.M, GM.C, GM.K, Wst, h)\(Fl*famp);

    GM.Rsc = max(abs(U0), min(abs(U0(U0~=0)))); 
%     UCs = PRECOCONT(@(Uw) GM.HBRESFUN(Uw, Fl*famp, h, Nt), U0, Wst, Wen, dw, Copt);
    UCs = CONTINUE(@(Uw) GM.HBRESFUN(Uw, Fl*famp, h, Nt), U0, Wst, Wen, dw, Copt);

    %% Plot
    Urms = sqrt([1, 0.5*ones(1,Nhc-1)]*UCs(1:end-1,:).^2);  % RMS Amplitude from Parseval's Thm.
    Uamp1 = abs(UCs(2,:)-1j*UCs(3,:));  % First Harmonic Amplitude
    Uph1 = rad2deg(angle(UCs(2,:)-1j*UCs(3,:)));  % First Harmonic Phase

    figure(1)
    subplot(2,2,1)
    plot(UCs(end,:), Uamp1, 'o-'); hold on
    ylabel('Displacement')

    subplot(2,2,2)
    semilogx(Urms, UCs(end,:), 'o-'); hold on 
    
    subplot(2,2,3)
    plot(UCs(end,:), Uph1, 'o-'); hold on
    ylabel('Phase (degs)')
    xlabel('Frequency (rad/s)')
end

%% EPMC Solution
As = -1;
Ae = 1.2;
da = 0.01;

% Copt.solverchoice = 2; 
% Copt.arclengthparm = 'K0NormalizedArcLength';


Uwx0 = zeros(Nhc, 1);
Uwx0(2:3) = 10^As/sqrt(2);  % "Mode" Initialization
Uwx0 = [Uwx0; sqrt(k/m); c/mul];
% [UwxC, dUwxC] = PRECOCONT(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt), Uwx0, As, Ae, da, Copt);
[UwxC, dUwxC] = CONTINUE(@(uwxa) GM.EPMCRESFUN(uwxa, Fl, h, Nt), Uwx0, As, Ae, da, Copt);

%% Plot EPMC Solution Over Freps
Ubb = (10.^UwxC(end, :)).*UwxC(1:end-3,:);
Urms = sqrt([1, 0.5*ones(1,Nhc-1)]*Ubb.^2);  % RMS Amplitude from Parseval's Thm.
Uamp1 = abs(Ubb(2,:)-1j*Ubb(3,:));  % First Harmonic Amplitude
Uph1 = rad2deg(angle(Ubb(2,:)-1j*Ubb(3,:)));  % First Harmonic Phase
    
figure(1)
subplot(2,2,1)
hold on;
plot(UwxC(end-2,:), Uamp1, 'k-', 'LineWidth', 1.2)
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