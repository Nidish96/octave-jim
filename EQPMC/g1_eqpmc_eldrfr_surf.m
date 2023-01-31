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
K = [2 -1;-1 2];
C = 0.01*K;

MDL = MDOFGEN(M, K, C, eye(2));

kt = 4;
muN = 2;
MDL = MDL.SETNLFUN(2+3, [1 0], @(t, u, varargin) JENKNL(t, u, kt, muN, varargin{:}), [], 4);
K0 = K+kt*[1 0;0 0]; % Linearized stiffness
[V, Wr] = eig(K0, M);
[Wr, si] = sort(sqrt(diag(Wr)));
V = V(:,si);

Nc = 2;  % Number of components
Nhmax = 3;  % Number of harmonics
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
Nt = 64;

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
Uwx0 = NSOLVE(@(X) MDL.EQPMCRESFUN([X; Xv(end)],  [1; 1], Fls, h, Nt, eps), Xv(1:end-1), opt);

% fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% Uwx = fsolve(@(X) MDL.EQPMCRESFUN([X; Xv(end)],  [1; 1], Fls, h, Nt, eps), Xv(1:end-1), fopt);

%% Do only ti = 1
Copt = struct('Nmax', 200, 'Display', 1, 'DynDscale', 0, 'solverchoice', 2, 'angopt', 1e2);
% Copt.Dscale = [kron([1e-2; 1e-1*ones(Nhc-1,1)], ones(2,1)); Wr; 1e-1*ones(2,1); 1.0];
% Copt.Dscale = [kron([1e-4; ones(Nhc-1,1)], ones(2,1)); 1e-2*ones(4,1); 1];

Astart = -2;
Aend = 1.35;  % 2
da = 0.1;

theta = pi/2/8;

% UwxLs = CONTINUE(@(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [cos(theta); sin(theta)], Fls, h, Nt, eps), ...
%     Uwx0, Astart, Aend, da, Copt);

% Sopt = struct('stepmax', 100, 'dynamicDscale', 0);
% UwxLs = solve_and_continue(Uwx0,@(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [cos(theta); sin(theta)], Fls, h, Nt, eps,1),...
%     Astart, Aend, da,Sopt);

% figure(2000)
% clf()
% % plot3((10.^UwxL{ti}(end,:))*cos(thetas(ti)), (10.^UwxL{ti}(end,:))*sin(thetas(ti)), UwxL{ti}(end-5+i,:), '.-', 'LineWidth', 2); hold on
% semilogx((10.^UwxLs(end,:))*cos(theta), UwxLs(end-4,:), '.-', 'LineWidth', 1); hold on
% grid on
% xlabel('Modal Amp $Q_1$')
% zlabel('Mode 1 Freq')

%% Continuation
% Copt = struct('Nmax', 200, 'Display', 1, 'DynDscale', 0, 'solverchoice', 2, 'angopt', 5e-2);
% Copt.Dscale = [kron([1e-2; 1e-1*ones(Nhc-1,1)], ones(2,1)); Wr; 1e-1*ones(2,1); 1.0];
% Copt.Dscale = [kron([1e-4; ones(Nhc-1,1)], ones(2,1)); 1e-2*ones(4,1); 1];

Copt = struct('Nmax', 500, 'Display', 1, 'DynDscale', 0);
% Copt.Dscale = [kron([1e-2; 1e-1*ones(Nhc-1,1)], ones(2,1)); Wr; 1e-1*ones(2,1); 1.0];
Astart = -2;
Aend = 2;  % 2
da = 0.05;
Copt.dsmax = 0.5;
Copt.dsmin = 0.05;

Sopt = struct('jac', 'full', 'stepmax', 2e3);
da = 0.05;

Nts = 32;
if analyze
    thetas = linspace(0, pi/2, Nts+2); thetas = thetas(2:end-1);
    UwxL = cell(size(thetas));
    for ti=1:length(thetas)
        try
%             UwxL{ti} = CONTINUE(@(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [cos(thetas(ti)); sin(thetas(ti))], Fls, h, Nt, eps), ...
%                 Uwx0, Astart, Aend, da, Copt);
            UwxL{ti} = solve_and_continue(Uwx0, @(Uwxl) MDL.EQPMCRESFUN(Uwxl,  [cos(thetas(ti)); sin(thetas(ti))], Fls, h, Nt, eps), ...
                Astart, Aend, da, Sopt);
        catch
            continue;
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
            '.-', 'LineWidth', 2); hold on
        catch
            continue;
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log')
    grid on
    xlabel('Modal Amplitude $Q_1$')
    ylabel('Modal Amplitude $Q_2$')
    xlim(10.^[Astart Aend])
    ylim(10.^[Astart Aend])
    set(gca, 'View', [-10 30])
    zlabel(sprintf('Mode %d Natural Frequency (rad/s)', i))

    if plotout
        set(gcf, 'Color', 'white')
        export_fig(sprintf('./FIGS/H%dW%d.png', Nhmax, i), '-dpng')
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
    if i==1
        zlim([0.95 1.45])
    else
        zlim([1 6])
    end

    if plotout
        set(gcf, 'Color', 'white')
        export_fig(sprintf('./FIGS/SURF_W%d.png',i), '-dpng')
    end
end

%% Compare against analytical expressions
falph = 0.5;
for i=1:2
    figure(100+i-1)
    clf()
    surf(Q1s, Q2s, WXs(:, :, i), 'EdgeColor', 'none','FaceAlpha',falph, 'FaceColor', 'b'); hold on
    surf(Q1s, Q2s, Wr(i)+3*0.5/8*((Q1s/2).^2+2*(Q2s/2).^2), 'EdgeColor', 'none', 'FaceAlpha', falph, 'FaceColor', 'r')
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlabel('Mode 1 Amp')
    ylabel('Mode 2 Amp')
%     xlim([min(Q1s(:)) max(Q1s(:))])
%     ylim([min(Q2s(:)) max(Q2s(:))])
    colormap(jet);
    set(gca, 'XTick', [1e-2 1e0], 'YTick', [1e-2 1e0])
    zlabel(sprintf('Mode %d Frequency (rad/s)', i))
    if i==1
        zlim([0.95 1.45])
    else
        zlim([1 6])
    end
    set(gca, 'View', [30 35])
    legend('Numerical', 'Analytical', 'Location', 'west')

    if plotout
        set(gcf, 'Color', 'white')
        export_fig(sprintf('./FIGS/SURF_W%d_comp.png',i), '-dpng')
    end    
end

%%
for j=1:2
    figure(12+j-1)
    clf()
%     surf(Q1s, Q2s, real(Zs(:, :, j, j)), imag(Zs(:, :, j, j)), 'EdgeColor', 'none'); hold on
%     colormap(jet)
%     yy = colorbar;
%     ylabel(yy, sprintf('Damping Factor $\\Im{\\zeta_{%d,%d}}$', j, j), 'Interpreter','latex')
    surf(Q1s, Q2s, real(Zs(:, :, j, j)), 'EdgeColor', 'none', 'FaceAlpha',falph); hold on
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlabel('Mode 1 Amp')
    ylabel('Mode 2 Amp')
    xlim([min(Q1s(:)) max(Q1s(:))])
    ylim([min(Q2s(:)) max(Q2s(:))])
    set(gca, 'XTick', [1e-2 1e0], 'YTick', [1e-2 1e0])
    zlabel(sprintf('Damping Coefficient ${c_{%d,%d}}$', j, j))
    if j==1
        set(gca, 'View', [-60 25])
    end
    if plotout
        set(gcf, 'Color', 'white')
        export_fig(sprintf('./FIGS/SURF_Z%d%d.png',j,j), '-dpng')
    end
end

figure(14)
clf()
% surf(Q1s, Q2s, real(Zs(:, :, 1, 2)), imag(Zs(:, :, 1, 2)), 'EdgeColor', 'none'); hold on
% colormap(jet);
% yy = colorbar;
% ylabel(yy, sprintf('Damping Factor $\\Im{\\zeta_{%d,%d}}$', 1, 2), 'Interpreter','latex')
surf(Q1s, Q2s, real(Zs(:, :, 1, 2)), 'EdgeColor', 'none','FaceAlpha',falph); hold on
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Mode 1 Amp')
ylabel('Mode 2 Amp')
xlim([min(Q1s(:)) max(Q1s(:))])
ylim([min(Q2s(:)) max(Q2s(:))])
set(gca, 'XTick', [1e-2 1e0], 'YTick', [1e-2 1e0])
zlabel(sprintf('Damping Coefficient $c_{%d,%d}$', 1, 2))
set(gca, 'View', [95 60])

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_Z%d%d.png',1,2), '-dpng')
end

figure(15)
clf()
% surf(Q1s, Q2s, real(Zs(:, :, 2, 1)), imag(Zs(:, :, 2, 1)), 'EdgeColor', 'none'); hold on
% colormap(jet);
% yy = colorbar;
% ylabel(yy, sprintf('Damping Factor $\\Im{\\zeta_{%d,%d}}$', 1, 2), 'Interpreter','latex')
surf(Q1s, Q2s, real(Zs(:, :, 2, 1)), 'EdgeColor', 'none','FaceAlpha',falph); hold on
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Mode 1 Amp')
ylabel('Mode 2 Amp')
xlim([min(Q1s(:)) max(Q1s(:))])
ylim([min(Q2s(:)) max(Q2s(:))])
set(gca, 'XTick', [1e-2 1e0], 'YTick', [1e-2 1e0])
zlabel(sprintf('Damping Coefficient $c_{%d,%d}$', 1, 2))
set(gca, 'View', [95 60])

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_Z%d%d.png',2,1), '-dpng')
end

%% Single Value (Contour) Plots
pint = 15;
colos = DISTINGUISHABLE_COLORS(Nlevs);

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
    wan = Wr(1)+3*0.5/8*((q1/2).^2+2*(q2/2).^2);
    semilogx(q1, wan, '-.', 'LineWidth', 2, 'Color', colos(i,:)); hold on
end
legend(aa, arrayfun(@(l) sprintf('$q_2=10^{%.2f}$', l), log10(levs(1:length(li))), 'UniformOutput', false), 'Location', 'northwest')
xlim(10.^[Astart Aend]);
xlabel('Mode 1 Amplitude')
ylabel('Mode 1 Frequency (rad/s)')
ylim([0.95 1.45])

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_WA%d%d.png',1,1), '-dpng')
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
    wan = Wr(1)+3*0.5/8*((q1/2).^2+2*(q2/2).^2);
    semilogx(q2, wan, '-.', 'LineWidth', 2, 'Color', colos(i,:)); hold on
end
legend(aa, arrayfun(@(l) sprintf('$q_1=10^{%.2f}$', l), log10(levs(1:length(li))), 'UniformOutput', false), 'Location', 'northwest')
xlim(10.^[Astart Aend]);
xlabel('Mode 2 Amplitude')
ylabel('Mode 1 Frequency (rad/s)')
ylim([0.95 1.45])

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_WA%d%d.png',2,1), '-dpng')
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
    wan = Wr(2)+3*0.5/8*((q1/2).^2+2*(q2/2).^2);
    semilogx(q1, wan, '-.', 'LineWidth', 2, 'Color', colos(i,:)); hold on
end
legend(aa, arrayfun(@(l) sprintf('$q_2=10^{%.2f}$', l), log10(levs(1:length(li))), 'UniformOutput', false), 'Location', 'northwest')
xlim(10.^[Astart Aend]);
xlabel('Mode 1 Amplitude')
ylabel('Mode 2 Frequency (rad/s)')
ylim([1 6])

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
    wan = Wr(2)+3*0.5/8*((q1/2).^2+2*(q2/2).^2);
    semilogx(q2, wan, '-.', 'LineWidth', 2, 'Color', colos(i,:)); hold on        
end
legend(aa, arrayfun(@(l) sprintf('$q_1=10^{%.2f}$', l), log10(levs(1:length(li))), 'UniformOutput', false), 'Location', 'northwest')
xlim(10.^[Astart Aend]);
xlabel('Mode 2 Amplitude')
ylabel('Mode 2 Frequency (rad/s)')
ylim([1 6])

if plotout
    set(gcf, 'Color', 'white')
    export_fig(sprintf('./FIGS/SURF_WA%d%d_eldrfr.png',2,2), '-dpng')
end

%%
ti = 1;

figure(2000)
clf()
% plot3((10.^UwxL{ti}(end,:))*cos(thetas(ti)), (10.^UwxL{ti}(end,:))*sin(thetas(ti)), UwxL{ti}(end-5+i,:), '.-', 'LineWidth', 2); hold on
plot((10.^UwxL{ti}(end,:))*cos(thetas(ti)), UwxL{ti}(end-5+i,:), '.-', 'LineWidth', 2); hold on
set(gca, 'XScale', 'log')
grid on
xlabel('Modal Amp $Q_1$')
zlabel('Mode 1 Freq')
