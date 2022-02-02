clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')

%%
m = 1;
c = 0.5;
k = 4;

Nc = 3;  % Number of components
Nhmax = 3;  % Number of harmonics
%% Harmonic Selection
hall = cell(1, Nc);
[hall{:}] = ndgrid(-Nhmax:Nhmax);
hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
h = hall(abs(sum(hall,2))<=Nhmax & sum(hall,2)>=0,:);

h(sum(h,2)==0 & h(:,1)<=0, :) = [];
h = [zeros(1,Nc); h];

figure(1)
clf()
plot(hall(:,1), hall(:,2), 'ko', 'MarkerFaceColor', 'w'); hold on
plot(h(:,1), h(:,2), '*'); hold on
grid on
axis equal
legend('All Harmonics', 'Selected Harmonic', 'Location', 'northoutside')

%% Setup Model
GM = MDOFGEN(m, k, c, 1.0);

%% Forcing
% ws = [pi sqrt(2)];
ws = pi.^(1:Nc);
rng(1);
hid = randi(length(h)-1,2,1)
hfrc = h(1+hid, :)
amps = [20; 30];
fext = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

%% Transient simulation
T0 = 0;
T1 = 30;

fsamp = 60;

opts = struct('Display', 'waitbar');
[T, U, Ud, Udd] = GM.HHTAMARCH(T0, T1, 1/fsamp, 0, 0, fext, opts);
% [T, U, Ud, Udd] = GM.HHTAMARCH(T0, T1, 1/fsamp, Ns(1,:)*X, Ns(1,:)*Xd, fext, opts);

%%
% sys = ss([0 1;-k/m -c/m], [0;1/m], [1 0], 0);
% [Y,T] = lsim(sys, fext(T), T);

%% MsHBM
Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));

E = QPHARMONICSTIFFNESS(GM.M, GM.C, GM.K, ws, h);
Fl = zeros(Nhc, 1);
Fl(1+(hid-1)*2+1) = amps;

D1 = QPHARMONICSTIFFNESS(0, 1, 0, ws, h);  % Time derivative matrix

X = E\Fl;
% X(2:3) = HARMONICSTIFFNESS(GM.M, GM.C, GM.K, hfrc*ws(:), 1)\Fl(2:3)
Xd = D1*X;
Xdd = D1^2*X;
Ns = QPTIMEINTERP(T(:).*ws, h);
%%
figure(2)
clf()
plot(T, U, 'o-')
grid on; hold on
% plot(T, Y, '.-');
plot(T, Ns*X, '-', 'LineWidth', 2)

figure(3)
clf()
plot(U, Ud, '.-'); grid on; hold on
plot(Ns*X, Ns*Xd, '-', 'LineWidth', 1)

%% Torus
if Nc==2
    Nt = 64;
    x   = reshape(QPTIMETRANS(X, h, Nt), repmat(Nt, 1, Nc));
    xd  = reshape(QPTIMETRANS(D1*X, h, Nt), repmat(Nt, 1, Nc));
    xdd = reshape(QPTIMETRANS(D1^2*X, h, Nt), repmat(Nt, 1, Nc));
    figure(4)
    clf()
    % plot3(x, xd, xdd, '-', 'LineWidth', 1)
    surf(x, xd, xdd, 'EdgeColor', 'none'); hold on
    plot3(U, Ud, Udd, 'k-')
end