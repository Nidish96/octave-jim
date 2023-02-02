clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')
addpath('../ROUTINES/export_fig/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

analyze = true;
plotout = true;
%% Parameters
M = eye(2);
K = [1 -1;-1 2];
W = sqrt(eig(K,M));

kt = 1;  % 4
muN = 1e-1;  % 2
K0 = K+kt*[1 0;0 0]; % Linearized stiffness
% C = 5e-3*K0;  % 0.01*K
C = 1e-2*M;

MDL = MDOFGEN(M, K, C, eye(2));
MDL = MDL.SETNLFUN(2+3, [1 0], @(t, u, varargin) JENKNL(t, u, kt, muN, varargin{:}), [], 4);
[V, Wr] = eig(K0, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

zts = diag(V'*C*V)./(2*Wr);
disp(zts)
disp([W Wr])

Nc = 2;  % Number of components
Nhmax = 1;  % Number of harmonics
%% Harmonic Selection
h = HSEL(Nhmax, [1 2]);
h = h(2:end, :);
Nhc = sum(all(h==0,2) + 2*any(h~=0,2));

%% Forcing
ws = [0.8; pi];
% rng(1);
% hid = randi(length(h)-1,2,1)
% hid = [1; 2];

hfrc = eye(2);

amps = 1.0*ones(Nc,1);  % 20
fext = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

[rinds0,zinds,hinds,rinds,iinds] = HINDS(2, h);

%% QP HB Simulation
Nt = 32;

[E, dEdw] = QPHARMONICSTIFFNESS(MDL.M, MDL.C, MDL.K, ws, h);
Fl = zeros(Nhc*2, 1);
Fl(rinds([1 3])) = 1.0;

X0 = E\Fl;

% opt = struct('Display', true, 'ITMAX', 100, 'crit', 30);
% X = NSOLVE(@(U) MDL.QPHBRESFUN([U; ws(:)], Fl, h, Nt, eps), X0, opt);

% fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% X = fsolve(@(U) MDL.QPHBRESFUN([U; ws(:)], Fl, h, Nt, eps), X0, fopt);
%% EQPMC Simulations
Fls = zeros(Nhc*2, 2);
Fls(rinds(1), 1) = 1;
Fls(rinds(3), 2) = 1;

xis = [M(:) K(:)]\C(:);  % Fit Proportional Constants (2 term Caughy series)

X0 = zeros(Nhc*2, 1);
X0([rinds(1:2) iinds(1:2)]) = kron([0;1], V(:,1));
X0([rinds(3:4) iinds(3:4)]) = kron([0;1], V(:,2));
Xv = [X0; Wr; xis; -1];
% [R, dRdU, dRdlA] = MDL.EQPMCRESFUN(Xv,  [1; 1], Fls, h, Nt, eps); 
% Uwx0 = NSOLVE(@(X) MDL.EQPMCRESFUN([X; Xv(end)],  [1; 1], Fls, h, Nt, eps), Xv(1:end-1), opt);

fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
Uwx0 = fsolve(@(X) MDL.EQPMCRESFUN([X; Xv(end)],  [1; 1], Fls, h, Nt, eps), Xv(1:end-1), fopt);

%% Continuation
% Copt = struct('Nmax', 200, 'Display', 1, 'DynDscale', 0, 'solverchoice', 2, 'angopt', 5e-2);
% Copt.Dscale = [kron([1e-2; 1e-1*ones(Nhc-1,1)], ones(2,1)); Wr; 1e-1*ones(2,1); 1.0];
% Copt.Dscale = [kron([1e-4; ones(Nhc-1,1)], ones(2,1)); 1e-2*ones(4,1); 1];

Copt = struct('Nmax', 100, 'Display', 1, 'DynDscale', 0, 'solverchoice', 1);
% Copt.Dscale = [kron([1e-2; 1e-1*ones(Nhc-1,1)], ones(2,1)); Wr; 1e-1*ones(2,1); 1.0];
Astart = -2;
Aend = 2;  % 2
da = 0.05;
Copt.dsmax = 0.25;
Copt.dsmin = 0.01;

Sopt = struct('jac', 'full', 'stepmax', 2e3, 'dynamicDscale', 1);
ds = 0.25;

Nts = 32;
if analyze
    thetas = linspace(0, pi/2, Nts+2); thetas = thetas(2:end-1);
    UwxL = cell(size(thetas));
    for ti=1:length(thetas)
        try
            UwxL{ti} = CONTINUE(@(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [cos(thetas(ti)); sin(thetas(ti))], Fls, h, Nt, eps), ...
                Uwx0, Astart, Aend, da, Copt);
%         UwxL{ti} = solve_and_continue(Uwx0, @(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [cos(thetas(ti)); sin(thetas(ti))], Fls, h, Nt, eps), ...
%             Astart, Aend, ds, Sopt);
        catch
            keyboard
        end

        if isempty(UwxL{ti}) || UwxL{ti}(end)<Aend
            sch = Copt.solverchoice;
            Copt.solverchoice = 3;
            UwxL{ti} = CONTINUE(@(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [cos(thetas(ti)); sin(thetas(ti))], Fls, h, Nt, eps), ...
                Uwx0, Astart, Aend, da, Copt);
            Copt.solverchoice = sch;
        end

        fprintf('Done %d/%d\n', ti, length(thetas));
    end
    save(sprintf('DATA/EQPSURFH%d_eldrfr.mat', Nhmax), 'thetas', 'UwxL')
else
    load(sprintf('DATA/EQPSURFH%d_eldrfr.mat', Nhmax), 'thetas', 'UwxL')
end

%% Plot 2D Surface
for i=1:2
    figure(1000+i-1)
    clf()
    for ti=1:length(thetas)
        try
        plot3((10.^UwxL{ti}(end,:))*cos(thetas(ti)), ...
            (10.^UwxL{ti}(end,:))*sin(thetas(ti)), UwxL{ti}(end-5+i,:), ...
            'o-', 'LineWidth', 2); hold on
        catch
            continue;
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log')
    grid on
    xlabel('Modal Amplitude $Q_1$')
    ylabel('Modal Amplitude $Q_2$')
%     xlim(10.^[Astart Aend])
%     ylim(10.^[Astart Aend])
    set(gca, 'View', [-10 30])
    zlabel(sprintf('Mode %d Natural Frequency (rad/s)', i))

    if plotout
        set(gcf, 'Color', 'white')
        export_fig(sprintf('./FIGS/H%dW%d_eldrfr.png', Nhmax, i), '-dpng')
    end
end
if Nhmax~=1
    return
end

%% Interpolate Surface
Npts = 50;
Q1s = zeros(Nts, Npts);
Q2s = zeros(Nts, Npts);
WXs = zeros(Nts, Npts, 4);
PHIs = zeros(Nts, Npts, MDL.Ndofs, 2);
Zs = zeros(Nts, Npts, 2, 2);

WXAs = cellfun(@(c) c(end-4:end,:), UwxL, 'UniformOutput', false);
P1s = cellfun(@(c) [c(MDL.Ndofs+(1:MDL.Ndofs),:)-1j*c(2*MDL.Ndofs+(1:MDL.Ndofs),:);c(end,:)], UwxL, 'UniformOutput', false);
P2s = cellfun(@(c) [c(3*MDL.Ndofs+(1:MDL.Ndofs),:)-1j*c(4*MDL.Ndofs+(1:MDL.Ndofs),:);c(end,:)], UwxL, 'UniformOutput', false);
for ti=1:length(thetas)
    WXAs{ti}(:, end) = interp1(WXAs{ti}(end,:), WXAs{ti}', Aend, 'pchip')';
    P1s{ti}(:, end) = interp1(P1s{ti}(end,:), P1s{ti}', Aend, 'pchip')';
    P2s{ti}(:, end) = interp1(P2s{ti}(end,:), P2s{ti}', Aend, 'pchip')';

    Q1s(ti, :) = logspace(Astart, Aend, Npts)*cos(thetas(ti));
    Q2s(ti, :) = logspace(Astart, Aend, Npts)*sin(thetas(ti));
    % Frequencies and Damping Coefficients
    WXs(ti, :, :) = interp1(WXAs{ti}(end,:), WXAs{ti}(1:end-1,:)', linspace(Astart, Aend, Npts));
    % Mode Shapes
    PHIs(ti, :, :, 1) = interp1(P1s{ti}(end,:), P1s{ti}(1:end-1,:)', linspace(Astart, Aend, Npts));  % Mode 1
    PHIs(ti, :, :, 2) = interp1(P2s{ti}(end,:), P2s{ti}(1:end-1,:)', linspace(Astart, Aend, Npts));  % Mode 2
    % Modal Damping Factors
    for i=1:2
        for j=1:2
            Zs(ti, :, i, j) = (WXs(ti, :, 3).'.*diag(conj(squeeze(PHIs(ti, :, :, i)))*M*squeeze(PHIs(ti, :, :, j)).') + ...
                WXs(ti, :, 4).'.*diag(conj(squeeze(PHIs(ti, :, :, i)))*K*squeeze(PHIs(ti, :, :, j)).')).*(WXs(ti, :, j)).';
        end
    end
end

%% Calculate Contours
Nlevs = 10;
levs = 10.^linspace(Astart, Aend, Nlevs+1);  levs=levs(1:end-1);

Conts = repmat(struct('Clevs', [], 'li', [], 'nis', []), 2,4);
for i=1:4
    Conts(1,i).Clevs = contour(Q1s, WXs(:,:,i), Q2s, levs);
    Conts(2,i).Clevs = contour(Q2s, WXs(:,:,i), Q1s, levs);
end

for qi=1:2
    for mj=1:4
        [~, Conts(qi,mj).li] = find(Conts(qi,mj).Clevs(1,:)==levs(:));
        Conts(qi,mj).nis = Conts(qi,mj).Clevs(2,Conts(qi,mj).li);
    end
end

%% Surface Plots
falph = 0.8;
for i=1:2
    figure(10+i-1)
    clf()
    surf(Q1s, Q2s, WXs(:, :, i), 'EdgeColor', 'none','FaceAlpha',falph); hold on
    for qi=1:2
        li = Conts(qi,i).li;
        nis = Conts(qi,i).nis;
        for k=1:length(li)
            if qi==1
                plot3(Conts(qi,i).Clevs(1, li(k)+(1:nis(k))), levs(k)*ones(nis(k),1), Conts(qi,i).Clevs(2, li(k)+(1:nis(k))), 'k-', 'LineWidth', 2);
            else
                plot3(levs(k)*ones(nis(k),1), Conts(qi,i).Clevs(1, li(k)+(1:nis(k))), Conts(qi,i).Clevs(2, li(k)+(1:nis(k))), 'w-', 'LineWidth', 2);
            end
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlabel('Mode 1 Amp')
    ylabel('Mode 2 Amp')
    xlim([min(Q1s(:)) max(Q1s(:))])
    ylim([min(Q2s(:)) max(Q2s(:))])
    colormap(jet);
    set(gca, 'XTick', [1e-2 1e0], 'YTick', [1e-2 1e0])
    zlabel(sprintf('Mode %d Frequency (rad/s)', i))
    set(gca, 'View', [120 20])

    if plotout
        set(gcf, 'Color', 'white')
        export_fig(sprintf('./FIGS/SURF_W%d_eldrfr.png',i), '-dpng')
    end
end

%% Single Value (Contour) Plots
pint = 15;
colos = DISTINGUISHABLE_COLORS(Nlevs);

load('./DATA/rdowndat_eldrfr.mat', 'tscale');
dpr = load('./DATA/rdowndat_eldrfr_mode2.mat', 'Amp', 'Freq');

figure(16)
clf()
% semilogx(Q1s(1:pint:end, :).', WXs(1:pint:end, :, 1).', '-', 'LineWidth', 2)
li = Conts(1,1).li;
nis = Conts(1,1).nis;
aa = gobjects(length(li),1);
for i=1:length(li)
    aa(i) = semilogx(Conts(1,1).Clevs(1,li(i)+(1:nis(i))), Conts(1,1).Clevs(2,li(i)+(1:nis(i))), '-', 'LineWidth', 2, 'Color', colos(i,:)); hold on

    q1 = Conts(1,1).Clevs(1,li(i)+(1:nis(i)));
    q2 = levs(i);
end
legend(aa, arrayfun(@(l) sprintf('$q_2=10^{%.2f}$', l), log10(levs(1:length(li))), 'UniformOutput', false), 'Location', 'northeast')
xlim(10.^[Astart Aend]);
xlabel('Mode 1 Amplitude')
ylabel('Mode 1 Frequency (rad/s)')
% ylim([0.95 1.45])

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_WA%d%d_eldrfr.png',1,1), '-dpng')
end

figure(17)
clf()
% semilogx(Q2s(1:pint:end, :).', WXs(1:pint:end, :, 1).', '-', 'LineWidth', 2)
li = Conts(2,1).li;
nis = Conts(2,1).nis;
for i=1:length(li)
    aa(i) = semilogx(Conts(2,1).Clevs(1,li(i)+(1:nis(i))), Conts(2,1).Clevs(2,li(i)+(1:nis(i))), '-', 'LineWidth', 2, 'Color', colos(i,:)); hold on

    q1 = levs(i);
    q2 = Conts(2,1).Clevs(1,li(i)+(1:nis(i)));
end
legend(aa, arrayfun(@(l) sprintf('$q_1=10^{%.2f}$', l), log10(levs(1:length(li))), 'UniformOutput', false), 'Location', 'northeast')
xlim(10.^[Astart Aend]);
xlabel('Mode 2 Amplitude')
ylabel('Mode 1 Frequency (rad/s)')
% ylim([0.95 1.45])

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_WA%d%d_eldrfr.png',2,1), '-dpng')
end

figure(18)
clf()
% semilogx(Q1s(1:pint:end, :).', WXs(1:pint:end, :, 2).', '-', 'LineWidth', 2)
li = Conts(1,2).li;
nis = Conts(1,2).nis;
for i=1:length(li)
    aa(i) = semilogx(Conts(1,2).Clevs(1,li(i)+(1:nis(i))), Conts(1,2).Clevs(2,li(i)+(1:nis(i))), '-', 'LineWidth', 2, 'Color', colos(i,:)); hold on

    q1 = Conts(1,2).Clevs(1,li(i)+(1:nis(i)));
    q2 = levs(i);
end
legend(aa, arrayfun(@(l) sprintf('$q_2=10^{%.2f}$', l), log10(levs(1:length(li))), 'UniformOutput', false), 'Location', 'northeast')
xlim(10.^[Astart Aend]);
xlabel('Mode 1 Amplitude')
ylabel('Mode 2 Frequency (rad/s)')
% ylim([1 6])

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_WA%d%d_eldrfr.png',1,2), '-dpng')
end

figure(19)
clf()
% semilogx(Q2s(1:pint:end, :).', WXs(1:pint:end, :, 2).', '-', 'LineWidth', 2)
li = Conts(2,2).li;
nis = Conts(2,2).nis;
for i=1:length(li)
    aa(i) = semilogx(Conts(2,2).Clevs(1,li(i)+(1:nis(i))), Conts(2,2).Clevs(2,li(i)+(1:nis(i))), '-', 'LineWidth', 2, 'Color', colos(i,:)); hold on

    q1 = levs(i);
    q2 = Conts(2,2).Clevs(1,li(i)+(1:nis(i)));
end
legend(aa, arrayfun(@(l) sprintf('$q_1=10^{%.2f}$', l), log10(levs(1:length(li))), 'UniformOutput', false), 'Location', 'northeast')
xlim(10.^[Astart Aend]);
xlabel('Mode 2 Amplitude')
ylabel('Mode 2 Frequency (rad/s)')

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_WA%d%d_eldrfr.png',2,2), '-dpng')
end

bb = plot(abs(dpr.Amp./(2*pi*dpr.Freq/tscale).^2/V(1,2)), (2*pi*dpr.Freq/tscale), 'ko-', ...
    'MarkerIndices', fix(logspace(0, log10(length(dpr.Amp)), 25)), ...
    'LineWidth', 2, 'MarkerFaceColor', 'w');
legend(bb, 'Transient');
legend([aa;bb], 'Location', 'northeast')

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_WA%d%d_eldrfr_wtr.png',2,2), '-dpng')
end
