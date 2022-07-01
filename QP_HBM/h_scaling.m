% clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

anim = false; 
plotfigs = true;
%%
m = 1;
% c = 0.5;
c = 0.01;
k = 1;

kt = 3;
muN = 0.5;

fnl = @(t, u, varargin) JENKNL(t, u, kt, muN, varargin{:});
% fnl = @(t,u,ud) deal(bt*u.^3, 3*bt*u.^2, zeros(size(u)));
% fnl = @(t,u,ud) deal(bt*u.^3 + c*ud, 3*bt*u.^2, c*ones(size(u)));

Nmtype = 3;  % Type of previous point jacobian construction

Nc = 2;  % Number of components
Nhmax = 7;  % Number of harmonics
%% Harmonic Selection
hall = cell(1, Nc);
[hall{:}] = ndgrid(-Nhmax:Nhmax);
hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
h = hall(sum(abs(hall),2)<=Nhmax & sum(hall,2)>=0,:);
% h = hall(abs(sum(hall,2))<=Nhmax & sum(hall,2)>=0,:);

h(sum(h,2)==0 & h(:,1)<=0, :) = [];
h = [zeros(1,Nc); h];

%% Setup Model
GM = MDOFGEN(m, k, c, 1.0);
GM = GM.SETNLFUN(2+3, 1.0, fnl, [], 1);

%% Forcing
% ws = [pi sqrt(2)];
% ws = [sqrt(k/m) pi.^(1:Nc-1)];
% ws = [pi.^(1:Nc)];
ws = sqrt(2:Nc+1);
% rng(1);
% hid = randi(length(h)-1,2,1)
% hid = [1; 2];

EYE = eye(Nc);
hid = zeros(Nc,1);
for i=1:Nc
    hid(i) = find(all(h==EYE(i,:),2))-1;
end
hfrc = h(1+hid, :);

amps = 5*ones(size(hid));  % forcing amplitude
fext = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

%% QP-HB solution
Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));

E = QPHARMONICSTIFFNESS(GM.M, GM.C, GM.K, ws, h);
Fl = zeros(Nhc, 1);
Fl(1+(hid-1)*2+1) = amps;

D1 = QPHARMONICSTIFFNESS(0, 1, 0, ws, h);  % Time derivative matrix

X0 = E\Fl;
Nt = 16;

%%
NTS = [16 32 64 128 256];
TTKS = zeros(length(NTS), 3);
opt = struct('Display', true, 'ITMAX', 200);
for qi=1:3
    GM.NLTs.qptype = qi;
    for ni=1:length(NTS)
        if qi<3 && ni>3
            TTKS(ni:end,qi) = nan;
            break
        end
            
        Nt = NTS(ni);
        
        tic
        for i=1:2
            X = NSOLVE(@(U) GM.QPHBRESFUN([U; 1], ws, Fl, h, Nt, eps), X0, opt);
        end
        TTKS(ni, qi) = toc/2;
    
        fprintf('============================\n')
        fprintf('%d: Done %d/%d\n', qi, ni, length(NTS))
        fprintf('============================\n')
    end
end

%%
figure(1)
clf()

colos = DISTINGUISHABLE_COLORS(3)*0.75;
aa = gobjects(3,1);
for ci=1:3
    aa(ci) = loglog(NTS, TTKS(:,ci), 's-', 'LineWidth', 2, ...
        'MarkerFaceColor', 'w', 'Color', colos(ci,:)); hold on
    legend(aa(ci), sprintf('Approach %d', ci))
end
legend(aa, 'Location', 'best')
xlabel('AFT Points (each dimension) $N_t$')
ylabel('CPU Time (s)')
set(gca, 'XTick', NTS)
xlim(NTS([1 end]))
grid on

set(gcf, 'Color', 'white')
if plotfigs
    export_fig('./FIGS/QPSCALING.png', '-dpng')
end