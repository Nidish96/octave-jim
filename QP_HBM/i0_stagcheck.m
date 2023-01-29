% QP marching approach for force estimation on staggered grid
clc
clear all
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

plotfigs = false;

Nmtype = 4;  % 4-full rotation; 5-stretched rotation
%% Frequency configuration
Nt = 32;
Nc = 2;
Nhmax = 5;

ws = sqrt(1:Nc);
% ws = [pi 1.0];

%% Harmonic Selection
hall = cell(1, Nc);
[hall{:}] = ndgrid(-Nhmax:Nhmax);
hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
h = hall(abs(sum(hall,2))<=Nhmax & sum(hall,2)>=0,:);

h(sum(h,2)==0 & h(:,1)<=0, :) = [];
h = [zeros(1,Nc); h];

%%
t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
deltau = t(2);
taus = cell(Nc, 1);
[taus{:}] = ndgrid(t);

tausn = cell(Nc, 1);
[tausn{:}] = ndgrid(1:Nt);

f = [0:Nt/2-1 -Nt/2:-1];
frqs = cell(Nc, 1);
[frqs{:}] = ndgrid(f);

%% Time-Transformation Matrix
switch Nc
    case 2
        A = 1/vecnorm(ws)*[ws(2) ws(1);-ws(1) ws(2)];
        A = [1 ws(1)/ws(2);0 1];
    otherwise 
        error('need to check');
end
Ai = inv(A);
ATi = inv(A');

%% 
% taup = A*tau
% w = A'*wp
taust = cell(Nc,1);
frqst = cell(Nc,1);
for ic=1:Nc
    taust{ic} = zeros(size(taus{ic}));
    frqst{ic} = zeros(size(frqs{ic}));
    for jc=1:Nc
        taust{ic} = taust{ic} + A(ic,jc)*taus{jc};
        frqst{ic} = frqst{ic} + A(ic,jc)*frqs{jc};
    end
end

%%
nor = [1; 2];
n = nor;
y = 0*2 + cos(n(1)*taus{1}+n(2)*taus{2});
yf = fft2(y);

%%
% nt = inv(A')*nor;
nt = nor;
yt = 0*2 + cos(nt(1)*taust{1}+nt(2)*taust{2});
ytf = fft2(yt);

%% Time Domain
figure(1)
clf()
surf(taus{1}, taus{2}, y, 'FaceColor', 'b'); hold on
surf(taust{1}, taust{2}, yt, 'FaceColor', 'r')

%% Frequency Domain
figure(2)
clf()
plot3(fftshift(frqs{1}), fftshift(frqs{2}), abs(fftshift(yf)), 'b.-'); hold on
plot3(fftshift(frqs{1}), fftshift(frqs{2}), abs(fftshift(ytf)), 'ro-')

box on; grid on;
% xlim([-1 1]*Nt/2)
% ylim([-1 1]*Nt/2)
xlabel('Frequency 1')
ylabel('Frequency 2')
zlabel('Component')