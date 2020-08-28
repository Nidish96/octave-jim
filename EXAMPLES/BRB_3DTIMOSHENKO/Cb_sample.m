clc
clear all
% addpath('../../ROUTINES/HARMONIC/')

%% Load Data
Nein = 8;
type = 'WGN';
DOF  = 'X';
fsamp = 2^18;

load(sprintf('./DATA/%dIN_%sRESP_%sDOFEX_samp%d.mat', Nein, type, DOF, log2(fsamp)), 'fsamp', ...
    'Ts', 'SensorLocs', 'ldof', 'Exc', 'Famps', 'Urecs', 'Udrecs', 'Uddrecs');
% DESCRIPTION:
%   fsamp       :  Sampling frequency used 2^18 Hz
%   Ts          :  3x1 Cell with time stamps: Ts{i} is Ntx1
%   SensorLocs  : 26x3 matrix with the locations of the 26 "sensors". Each row is [X, Y, Z]
%   ldof        : integer indicating the DOF that the forcing was applied to (1-sensor 1 x; 2-sensor 1 y; 3-sensor 1 z; 4-sensor 2 x; ...)
%   Exc         : 3x1 Cell of Ntx1 time series of the forcing that is applied at ldof (location can be obtained from SensorLocs)
%   Famps       : 3x1 vector of the maximum amplitude of the White Gaussian Noise Spectrum in the Frequency Domain
%   Urecs       : 3x1 Cell of 54xNt Displacement response at the 3 directions of the 26 sensors ([u1x;u1y;u1z;u2x;u2y;u2z;...])
%   Udrecs      : 3x1 Cell of 54xNt Velocity response of the 54 chosen DOFs
%   Uddrecs     : 3x1 Cell of 54xNt Acceleration response of the 54 chosen DOFs


%% Example time series plotting
sensor_idx = 1;  % Sensor ID
sensor_coords = SensorLocs(sensor_idx, :);  % [X Y Z] coordinates of chosen sensor

sol = Uddrecs;

XL = [0 3e-1];


r_idx = (sensor_idx-1)*3;  % Index pointing to location in the Urecs{i} data
% Ts = {Ts{1}};
figure(1)
clf()
set(gcf, 'Position', [2721 184 560 740])
aa = gobjects(size(Ts));
for fi=length(Ts):-1:1
    subplot(4,1,1)
    plot(Ts{fi}, sol{fi}(r_idx+1, :), '-'); hold on
    ylabel('AX')
    xlim(XL)
    
    subplot(4,1,2)
    plot(Ts{fi}, sol{fi}(r_idx+2, :), '-'); hold on
    ylabel('AY')
    xlim(XL)
    
    subplot(4,1,3)
    plot(Ts{fi}, sol{fi}(r_idx+3, :), '-'); hold on
    ylabel('AZ')
    xlim(XL)
    
    subplot(4,1,4)
    aa(fi) = plot(Ts{fi}, Exc{fi}, '-'); hold on
	legend(aa(fi), sprintf('F=%d N', Famps(fi)))
    ylabel('Force')
    xlim(XL)
end
subplot(4,1,4)
xlabel('Time (s)')
legend(aa(1:end), 'Location', 'best')

%% Frequency Domain response - 3 DOFs (x, y, z) of sensor 1
% No windowing is done, so the frequency domain looks pretty noisy
figure(2)
clf()
set(gcf, 'Position', [3281 183 560 740]);
aa = gobjects(size(Famps));
for fi=length(Famps):-1:1
    [freqs, UF] = FFTFUN(Ts{fi}, [Uddrecs{fi}(1:3, :)']);
    [~, FF] = FFTFUN(Ts{fi}, Exc{fi}');

    subplot(4, 1, 1)
    semilogy(freqs, abs(UF(:,1))); hold on
    xlim([0 3e3])
    ylabel('AX')
    
    subplot(4, 1, 2)
    semilogy(freqs, abs(UF(:,2))); hold on
    xlim([0 3e3])
    ylabel('AY')
    
    subplot(4, 1, 3)
    semilogy(freqs, abs(UF(:,3))); hold on
    xlim([0 3e3])
    ylabel('AZ')
    
    subplot(4, 1, 4)
    aa(fi) = semilogy(freqs, abs(FF)); hold on
    legend(aa(fi), sprintf('F=%d N', Famps(fi)))
    xlim([0 3e3])
    ylim([1e0 1e4])
    ylabel('Force')
end
subplot(4, 1, 4)
xlabel('Frequency (Hz)');
legend(aa(1:end), 'Location', 'best')

%% Plot Location of sensors
ti = 500;  % Choose frame to depict

f_idx = 1;

sc = 1e3;  % Scaling displacements to amplify response visually

figure(3)
clf()

% Plot Sensor Locations + displacements
plot3(SensorLocs(:, 1)+sc*Urecs{f_idx}(1:3:end,ti), ...
    SensorLocs(:, 2)+sc*Urecs{f_idx}(2:3:end,ti), ...
    SensorLocs(:, 3)+sc*Urecs{f_idx}(3:3:end,ti), ...
    'b.', 'MarkerSize', 10); hold on
% Beam Neutral Axis
Xna = SensorLocs([1; 2; 6; 10; 14; 18; 22; 26], 1);
Ux = 0.25*(Urecs{f_idx}([1; 2; 6; 10; 14; 18; 22; 26]*3-3+1,ti)+...
    Urecs{f_idx}([1; 3; 7; 11; 15; 19; 23; 26]*3-3+1,ti)+...
    Urecs{f_idx}([1; 4; 8; 12; 16; 20; 24; 26]*3-3+1,ti)+...
    Urecs{f_idx}([1; 2; 9; 13; 17; 21; 25; 26]*3-3+1,ti));
Uy = 0.25*(Urecs{f_idx}([1; 2; 6; 10; 14; 18; 22; 26]*3-3+2,ti)+...
    Urecs{f_idx}([1; 3; 7; 11; 15; 19; 23; 26]*3-3+2,ti)+...
    Urecs{f_idx}([1; 4; 8; 12; 16; 20; 24; 26]*3-3+2,ti)+...
    Urecs{f_idx}([1; 2; 9; 13; 17; 21; 25; 26]*3-3+2,ti));
Uz = 0.25*(Urecs{f_idx}([1; 2; 6; 10; 14; 18; 22; 26]*3-3+3,ti)+...
    Urecs{f_idx}([1; 3; 7; 11; 15; 19; 23; 26]*3-3+3,ti)+...
    Urecs{f_idx}([1; 4; 8; 12; 16; 20; 24; 26]*3-3+3,ti)+...
    Urecs{f_idx}([1; 2; 9; 13; 17; 21; 25; 26]*3-3+3,ti));
plot3(Xna+sc*Ux, sc*Uy, sc*Uz, 'ko-'); 
plot3(Xna, Xna*0, Xna*0, 'k--');

axis equal
grid on

xlim([-0.1 0.8])
% zlim(0.1*[-1 1])

xlabel('X')
ylabel('Y')
zlabel('Z')

title(sprintf('Frame %d. Blue points are sensors, black line is neutral axis. Joint not depicted.', ti))
pause(0.1)