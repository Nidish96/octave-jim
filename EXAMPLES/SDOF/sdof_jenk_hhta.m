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

kt = 5.0*mul;
muN = 0.1*mul;

fnl = @(t, u, varargin) JENKFORCE(t, u, kt, muN, varargin{:});

%% Setup Model
GM = MDOFGEN(m, k ,c, 1.0);
GM = GM.SETNLFUN(2+3, 1.0, fnl);  % See MDOFGEN class for understanding meaning of the arguments

%% Haversine impulse
pfq = sqrt(k/m)*2;
famp = 1.0;
fext = @(t) famp*sin(2*pi*pfq*t).^2.*(t<=1/pfq & t>=1/pfq/2);

%% Transient
T0 = 0;
T1 = 1/(c/2)*10;

fsamp = 30;

opts = struct('Display', 'waitbar');
[T, U, Ud, Udd] = GM.HHTAMARCH(T0, T1, 1/fsamp, 0, 0, fext, opts);

figure(10)
clf()
subplot(3,1,1)
plot(T, U, '.-');
ylabel('Displ.')
subplot(3,1,2)
plot(T, Ud, '.-');
ylabel('Vel.')
subplot(3,1,3)
plot(T, Udd, '.-');
ylabel('Accel.')
xlabel('Time (s)')