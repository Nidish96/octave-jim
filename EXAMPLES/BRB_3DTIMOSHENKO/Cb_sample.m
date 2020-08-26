clc
clear all

%% Load Data
Nein = 8;
type = 'WGN';
DOF  = 'Z';

load(sprintf('./DATA/%dIN_%sRESP_%sDOFEX.mat', Nein, type, DOF), 'fsamp', ...
    'Ts', 'SensorLocs', 'ldof', 'Exc', 'Famps', 'Urecs', 'Udrecs', 'Uddrecs');
% load(sprintf('./DATA/new/%dIN_%sRESP_%sDOFEX.mat', Nein, type, DOF), 'fsamp', ...
%     'Ts', 'SensorLocs', 'ldof', 'Exc', 'Famps', 'Urecs', 'Udrecs', 'Uddrecs');
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
DOF_idx    = 1;  % DOF ID (1-X, 2-Y, 3-Z)
sensor_coords = SensorLocs(sensor_idx, :);  % [X Y Z] coordinates of chosen sensor

r_idx = (sensor_idx-1)*3+DOF_idx;  % Index pointing to location in the Urecs{i} data
% Ts = {Ts{1}};
figure(1)
clf()
aa = gobjects(size(Ts));
for i=1:length(Ts)
    
    aa(i) = plot(Ts{i}, Urecs{i}(r_idx, :), '-'); hold on
%     aa(i) = plot(Ts{i}, Exc{i}(:), '-'); hold on
    legend(aa(i), sprintf('%d N', Famps(i))) 
end
legend(aa(1:end), 'Location', 'northeast')

xlabel('Time (s)')
ylabel('Displacement (m)')

return
%% Plot Location of sensors
timei = 500/fsamp;

ti = fix(timei*fsamp);
f_idx = 1;

sc = 1e3;  % Scaling displacements to amplify response visually

for ti = 1:1000:100000
figure(2)
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
zlim(0.1*[-1 1])

xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')

title(sprintf('Frame %d', ti))
pause(0.1)
end