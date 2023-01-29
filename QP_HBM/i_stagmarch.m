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
Nt = 16;
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

%% Time-Transformation Matrix
switch Nc
    case 2
        A = 1/vecnorm(ws)*[ws(2) ws(1);-ws(1) ws(2)];
    otherwise 
        error('need to check');
end
Ai = inv(A);

%%
taust = cell(Nc,1);
for ic=1:Nc
    taust{ic} = zeros(size(taus{ic}));
    for jc=1:Nc
        taust{ic} = taust{ic} + A(ic,jc)*taus{jc};
    end
end

%%
figure(1)
clf()
plot(taus{1}, taus{2}, 'k.'); hold on
plot(taust{1}, taust{2}, 'r.')

%% Displacement (first harmonic displacement for each component)
hid = eye(Nc);

amps = (1+(1:Nc));
utau = amps(1)*cos(taus{1});
for ti=2:Nc
    utau = utau + amps(ti)*cos(taus{ti});
end
ut = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

%% Elastic Dry Friction Parameters
kt = 1.0;
% muN = sqrt(sum(amps.^2)/2)/2;
switch Nc
    case 1
        muN = 1;
    otherwise
        muN = 2.5;
end

%%
