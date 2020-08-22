clc
clear all
addpath('../../ROUTINES/HARMONIC/')
addpath('../../ROUTINES/SOLVERS/')

%% Parameters
m = 1.00;   % mass
c = 0.20;   % linear damping
k = 4*pi^2; % linear stiffness
al = 5.0;   % cubic stiffness coefficient

%% MHBM FRF
Nt   = 256;                      % Number of time points per cycle
h    = [0 1 3];            % List of harmonics to balance
Nhc  = sum(h==0)+2*sum(h~=0);
Fl   = [0 1 0 zeros(1, Nhc-3)]';

Ws   = 2*pi*0.5;      % Starting Frequency (rad/s)
We   = 2*pi*1.5;  % Ending Frequency (rad/s)
ds   = 0.1;      % Arc-length step size

F_level = [0.5 1.0 2.0 4.0];
U = cell(size(F_level));

opts = struct('reletol', 1e-10, 'rtol', 1e-6, 'utol', 5e-6, 'etol', ...
                  1e-6, 'ITMAX', 20, 'Display', false, 'Dscale', ...
                  ones(Nhc,1));
Copt = struct('Nmax', 1000, 'Display', 1, 'angopt', 1e-2, 'opts', ...
                  opts);

for fi = 1:length(F_level)
    Elin = HARMONICSTIFFNESS(m, c, k, Ws, h);
    U0   = Elin\(Fl*F_level(fi));  % Initial Linear Prediction

%     Copt.Dscale = [ones(length(U0),1); (Ws+We)/2];

    U{fi} = CONTINUE(@(Uw) SDOF_NL_HBRESFUN(Uw, m, c, k, Fl*F_level(fi), ...
                               @(t,x,xd) DUFFNL(t,x,xd,al), h, Nt), ...
                U0, Ws, We, ds, Copt);
end

%% MHBM EPMC (Backbones through Extended periodic Motion Concept)
As = -1;
Ae = +1;
da = 0.01;

U0 = [0; 0; 1; zeros(Nhc-3,1); sqrt(k/m); c/m];

Copt.opts.Dscale = ones(Nhc+2,1);
UwxA = CONTINUE(@(Uwxa) SDOF_NL_EPMCRESFUN(Uwxa, m, c, k, Fl, @(t, x, xd) DUFFNL(t,x,xd,al), h, Nt), ...
    U0, As, Ae, da, Copt);


%% Plot FRF
figure(1)
clf()

figure(2)
clf()

for fi = 1:length(F_level)
    figure(1)
    Nc = size(U{fi}, 2);
    plot(U{fi}(end, :)/2/pi, sqrt(sum(U{fi}(1:end-1, :).^2.*repmat([1; 0.5*ones(Nhc-1,1)], 1, Nc), 1)), '-'); hold on
    
    figure(2)
    plot(U{fi}(end, :)/2/pi, rad2deg(angle(U{fi}(2, :)-1j*U{fi}(3,:))), '-'); hold on
end

figure(1)
Nc = size(UwxA,2);
plot(UwxA(end-2,:)/2/pi, sqrt(sum((UwxA(1:end-3,:)).^2.*repmat([1; 0.5*ones(Nhc-1,1)], 1, Nc), 1)).*10.^UwxA(end,:), 'k--');

xlim([Ws We]/2/pi)
xlabel('Frequency (Hz)')
ylabel('RMS Displacement Amplitude (m)')

figure(2)
xlabel('Frequency (Hz)')
ylabel('Response Phase (degs)')

figure(3) 
clf()
semilogx(sqrt(sum((UwxA(1:end-3,:)).^2.*repmat([0; ones(Nhc-1,1)], 1, Nc), 1)).*10.^UwxA(end,:), UwxA(end-2,:)/2/pi, 'k-')
xlabel('Displacement Amplitude (m)')
ylabel('EPMC Natural Frequency (Hz)')

figure(4) 
clf()
semilogx(sqrt(sum((UwxA(1:end-3,:)).^2.*repmat([0; ones(Nhc-1,1)], 1, Nc), 1)).*10.^UwxA(end,:), UwxA(end-1,:)./(2*UwxA(end-2,:)), 'k-')
xlabel('Displacement Amplitude (m)')
ylabel('EPMC Damping Factor')