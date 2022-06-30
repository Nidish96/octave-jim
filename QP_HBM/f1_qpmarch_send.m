clc
clear all

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

Nmtype = 2;  % "Approaches" for constructing N matrix : 1-interpolation ; 2-fdm
%% Frequency configuration
Nt = 16;  % Number of time points per time scale
Nc = 2;   % Number of frequency components
Nhmax = 5;% Max harmonic

ws = sqrt(1:Nc);  % Frequency Scales

%% Harmonic Selection
hall = cell(1, Nc);
[hall{:}] = ndgrid(-Nhmax:Nhmax);
hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
h = hall(abs(sum(hall,2))<=Nhmax & sum(hall,2)>=0,:);

h(sum(h,2)==0 & h(:,1)<=0, :) = [];
h = [zeros(1,Nc); h];

%% Construct Meshgrids for hyper-time domain
t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
deltau = t(2);
taus = cell(Nc, 1);
[taus{:}] = ndgrid(t);

tausn = cell(Nc, 1);
[tausn{:}] = ndgrid(1:Nt);

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
switch Nc
    case 1
        muN = 1;
    otherwise
        muN = 2.5;
end

%% MAIN PORTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Approach 1: Build N matrix with Lagrange Shape Functions
dt_vec = ws*deltau/vecnorm(ws);  % vector corresponding to deltatau amplitude in real time dxn

ptsb = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^Nc-1))));  % using binary for the construction of points on a unit square
oppi = (2^Nc):-1:1;  % diagonally opposite points are retrieved using the binary inverses
xis = ptsb*(-deltau);  % coordinates in "tau" space relative to origin

Lm = deltau^Nc;  % Lebesgue Measure of each cell in tau space
Nsf = prod(abs(xis(oppi,:)-(-dt_vec)), 2)'/Lm;  % Lagrange Shape Functions to interpolate previous point in the cell diagonally behind (in tau space)
    
ijs = cell2mat(cellfun(@(c) c(:), tausn, 'UniformOutput', false)');  % indices of all points
evijs = mod(repmat(ijs, 1, 1, 2^Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the cell that is diagonally behind
evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

% Build sparse interpolation matrix
tic
Nmat1 = sparse(repmat((1:Nt^Nc)', 1, 2^Nc), evns, repmat(Nsf, Nt^Nc, 1));
toc

%% Scheme 2: Build N matrix from FD discretization
ptsb = eye(Nc);

Nsf = ws/sum(ws);
    
ijs = cell2mat(cellfun(@(c) c(:), tausn, 'UniformOutput', false)');  % indices of all points
evijs = mod(repmat(ijs, 1, 1, Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the cell that is diagonally behind
evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

% Build sparse interpolation matrix
tic
Nmat2 = sparse(repmat((1:Nt^Nc)', 1, Nc), evns, repmat(Nsf, Nt^Nc, 1));
toc

%% Choose Nmat
switch Nmtype
    case 1
        Nmat = Nmat1;
    case 2
        Nmat = Nmat2;
end

%% Implicit "March"
cstick = @(f) abs(kt*(speye(Nt^Nc)-Nmat)*utau(:) + Nmat*f)<muN;  % stick case
cslip = @(f) abs(kt*(speye(Nt^Nc)-Nmat)*utau(:) + Nmat*f)>=muN;  % slip case

fspr = @(f) f - (kt*(speye(Nt^Nc)-Nmat)*utau(:)+Nmat*f).*cstick(f) - (muN*sign(kt*(speye(Nt^Nc)-Nmat)*utau(:)+Nmat*f)).*cslip(f);
fsp = @(f) deal(f - (kt*(speye(Nt^Nc)-Nmat)*utau(:)+Nmat*f).*cstick(f) - (muN*sign(kt*(speye(Nt^Nc)-Nmat)*utau(:) + Nmat*f)).*cslip(f), ...
    speye(Nt^Nc) - (Nmat).*cstick(f));

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
fsol = fsolve(fsp, kt*utau(:), opt);

% opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display', 'iter');
% fsol = fsolve(fspr, kt*utau(:), opt);
%% END OF MAIN PORTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Temporal verification using Backward Euler 
if Nc<3
    Tmax = 2*pi/min(ws)*100;
else
    Tmax = 2*pi/min(ws)*200;
end
dt = 0.001;

tvec = (0:dt:Tmax)';
us = cos(tvec.*ws)*amps(:);
ft = zeros(length(us),1);
for ti=2:length(us)
    fspred = kt*(us(ti)-us(ti-1)) + ft(ti-1);
    if abs(fspred)<muN
        ft(ti) = fspred;
    else
        ft(ti) = muN*sign(fspred);
    end
end

%% Reconstruct f from fsol (in hyper-time space)
fmx = reshape(fsol, [repmat(Nt, 1, Nc) ones(1,Nc==1)]);  % f in matrix form (in tau space)

repinds = mat2cell(repmat([1:Nt 1]', 1, Nc), Nt+1, ones(1, Nc));
fmxp = fmx(repinds{:});
tausp = cell(Nc, 1);
[tausp{:}] = ndgrid([t; 2*pi]);

Xq = mat2cell(wrapTo2Pi(tvec.*ws), length(tvec), ones(1, Nc));
fq = interpn(tausp{:}, fmxp, Xq{:}, 'cubic');
%%
figure((Nmtype-1)*2+1)
clf(); 
plot(tvec, ft, '-', 'LineWidth', 2); hold on
plot(tvec, fq, '-.', 'LineWidth', 2)
switch Nc
    case 1
        xlim([0 20])
    case 2
        xlim([0 30])
    case 3
        xlim([0 30])
end

ll = legend('Transient March', 'Hyper-Time March', 'Location', 'southoutside');
set(ll, 'NumColumns', 2)

xlabel('Time (s)')
ylabel('Force (N)')
set(gcf, 'Color', 'white')

%%
figure((Nmtype-1)*2+2)
clf()
set(gcf, 'Position', [400 500 900 450])

subplot(1,2,1)
switch Nc
    case 2
        scatter3(Xq{1}, Xq{2}, ft, ones(size(ft)), ft)
    case 3
        scatter3(Xq{:}, [], ft, '.')
end
title('Transient March')
colormap(jet);
xx = colorbar('northoutside');
xlabel(xx, 'Force (N)')

subplot(1,2,2)
switch Nc
    case 2
        surf(tausp{:}, fmxp)
    case 3
        scatter3(tausp{1}(:), tausp{2}(:), tausp{3}(:), [], fmxp(:), 'filled')
end
title('Hyper-Time March')
colormap(jet)
xx = colorbar('northoutside');
xlabel(xx, 'Force (N)')
for i=1:2
    subplot(1,2,i)
    tics = (0:.25:1);
    ticls = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
    set(gca, 'XTick', tics*2*pi)
    set(gca, 'XTickLabel', ticls)
    set(gca, 'YTick', tics*2*pi)
    set(gca, 'YTickLabel', ticls)
    xlabel('Time $\tau_1$', 'Rotation', 28)
    ylabel('Time $\tau_2$', 'Rotation', -38)
    switch Nc
        case 2
            zlabel('Force (N)')
        case 3
            zlabel('Time $\tau_3$')
            set(gca, 'ZTick', tics*2*pi)
            set(gca, 'ZTickLabel', ticls)
    end
end
set(gcf, 'Color', 'white')
