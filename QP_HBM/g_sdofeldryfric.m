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

Nc = 3;  % Number of components
Nhmax = 7;  % Number of harmonics
%% Harmonic Selection
hall = cell(1, Nc);
[hall{:}] = ndgrid(-Nhmax:Nhmax);
hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
h = hall(sum(abs(hall),2)<=Nhmax & sum(hall,2)>=0,:);
% h = hall(abs(sum(hall,2))<=Nhmax & sum(hall,2)>=0,:);

h(sum(h,2)==0 & h(:,1)<=0, :) = [];
h = [zeros(1,Nc); h];

switch Nc
    case 2
        figure(1)
        clf()
        plot(hall(:,1), hall(:,2), 'ko', 'MarkerFaceColor', 'w'); hold on
        plot(h(:,1), h(:,2), '*'); hold on
        grid on
        xlabel('Index $h_1$')
        ylabel('Index $h_2$')
        axis equal
        legend('All Harmonics', 'Selected Harmonic', 'Location', 'northoutside')
        set(gcf, 'Color', 'white')
        if plotfigs
            export_fig('./FIGS/G_selharms.png', '-dpng')
        end
end

%% Setup Model
GM = MDOFGEN(m, k, c, 1.0);
GM = GM.SETNLFUN(2+3, 1.0, fnl, [], Nmtype);

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

%% Transient simulation
T0 = 0;
T1 = 120;

fsamp = 60;

opts = struct('Display', 'waitbar');
[T, U, Ud, Udd] = GM.HHTAMARCH(T0, T1, 1/fsamp, 0, 0, fext, opts);

%% QP-HB solution
Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));

E = QPHARMONICSTIFFNESS(GM.M, GM.C, GM.K, ws, h);
Fl = zeros(Nhc, 1);
Fl(1+(hid-1)*2+1) = amps;

D1 = QPHARMONICSTIFFNESS(0, 1, 0, ws, h);  % Time derivative matrix

X0 = E\Fl;
Nt = 16;

% [R0, dR0] = GM.QPHBRESFUN([X0; 1], ws, Fl, h, Nt, eps);
% %%

tic
opt = struct('Display', true, 'ITMAX', 200);
X = NSOLVE(@(U) GM.QPHBRESFUN([U; 1], ws, Fl, h, Nt, eps), X0, opt);
toc

% fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% X = fsolve(@(U) GM.QPHBRESFUN([U; 1], ws, Fl, h, Nt, eps), X0, fopt);

Xd = D1*X;
Xdd = D1^2*X;
Ns = QPTIMEINTERP(T(:).*ws, h);

%%
figure(2)
clf()
set(gcf, 'Position', [400 500 900 450])
set(gcf, 'Color', 'white')

aa = plot(T, U, 'b-', 'LineWidth', 2);
grid on; hold on
bb = plot(T, Ns*X, 'r-', 'LineWidth', 2);

xlabel('Time (s)')
ylabel('Displacement (m)')
yl = ylim;
switch Nc
    case 1
        xl = [0 20];
    case 2
        xl = [0 20];
    case 3
        xl = [0 20];
end
mm=mesh([xl(:) xl(:)], [yl;yl], zeros(2), 'FaceColor', 'none', 'EdgeColor', 'k');
ylim(yl); 

ll = legend([aa bb], 'Transient', 'QP-HBM', 'Location', 'northoutside');
set(ll, 'NumColumns', 2);
if plotfigs
%     delete(mm)
    export_fig(sprintf('./FIGS/G_ut_%d_nm%d.png', Nc, Nmtype), '-dpng')

    xlim(xl)
    ylim(yl)
    export_fig(sprintf('./FIGS/G_utzoom_%d_nm%d.png', Nc, Nmtype), '-dpng')
end

%%
figure(3)
clf()
plot(U, Ud, 'b-', 'LineWidth', 2); grid on; hold on
plot(Ns*X, Ns*Xd, 'r-', 'LineWidth', 1)
ll = legend('Transient', 'QP-HBM', 'Location', 'northeast');

xlabel('Displacement ($m$)')
ylabel('Velocity ($m/s$)')
set(gcf, 'Color', 'white')
if plotfigs
    export_fig(sprintf('./FIGS/G_uv_%d_nm%d.png', Nc, Nmtype), '-dpng')
end
%% Torus
if Nc==2 || Nc==3
    Nt = 64;
    x   = reshape(QPTIMETRANS(X, h, Nt), repmat(Nt, 1, Nc));
    xd  = reshape(QPTIMETRANS(D1*X, h, Nt), repmat(Nt, 1, Nc));
    xdd = reshape(QPTIMETRANS(D1^2*X, h, Nt), repmat(Nt, 1, Nc));
    
%     x   = [x   x(:,1)  ; x(1,:)   x(1,1)];
%     xd  = [xd  xd(:,1) ; xd(1,:)  xd(1,1)];
%     xdd = [xdd xdd(:,1); xdd(1,:) xdd(1,1)];
    figure(4)
    clf()
    Nlast = 1;
    plot3(U(Nlast:end), Ud(Nlast:end), Udd(Nlast:end), 'b-', 'LineWidth', 2)
    grid on; hold on
    
    % plot3(x, xd, xdd, '-', 'LineWidth', 1)
    if Nc==2
        surf(x, xd, xdd, 'EdgeColor', 'k'); hold on
    elseif Nc==3
        scatter3(x(:), xd(:), xdd(:), [], xdd(:), '.')
    end
    
    colormap(jet);
    xx = colorbar('southoutside');
    xlabel(xx, 'Acceleration ($m/s^2$)', 'Interpreter', 'latex')
    colormap(jet)    
    
    xlabel('Displacement ($m$)', 'Rotation', 15)
    ylabel('Velocity ($m/s$)', 'Rotation', -24)
    zlabel('Acceleration ($m/s^2$)')
    
    ll = legend('Transient', 'QP-HBM', 'Location', 'northeast');
    set(gca, 'View', [-22 85])
    set(gcf, 'Color', 'white')
    if plotfigs
        export_fig(sprintf('./FIGS/G_torus_%d_nm%d.png', Nc, Nmtype), '-dpng')
    end
end

%% Frequency Content
% wndw = ones(length(T), 1);
wndw = blackmanharris(length(T));
% wndw = hanning(length(T));

[frqs, uf] = FFTFUN(T', U'.*wndw);
[~, xf] = FFTFUN(T', (Ns*X).*wndw);

figure(5)
clf()
set(gcf, 'Color', 'white')

semilogy(2*pi*frqs, abs(uf), 'b-'); hold on
% plot(frqs, abs(xf), 'r-')
stem(abs(h*ws(:)), abs([X(1); X(2:2:end)+1j*X(3:2:end)]), 'filled', 'k')
xlabel('Frequency (rad/s)')
ylabel('Displacement (m)')
yl = ylim;
switch Nc
    case 1
        xl = [0 2]*2*pi;
        yl = [1e-7 1e1];
    case 2
        xl = [0 1.5]*2*pi;
        yl = [1e-6 1e1];
    case 3
        xl = [0 10];
        yl = [1e-4 1e1];
end
mesh([xl(:) xl(:)], [yl;yl], zeros(2), 'FaceColor', 'none', 'EdgeColor', 'k')
legend('Transient', 'QP-HBM')
if plotfigs
    export_fig(sprintf('./FIGS/G_fcont_%d_nm%d.png', Nc, Nmtype), '-dpng')
    
    xlim(xl)
    ylim(yl)
    export_fig(sprintf('./FIGS/G_fcontzoom_%d_nm%d.png', Nc, Nmtype), '-dpng')
end

%% Animate
if Nc==2 && anim
    Nt = 64;
    x   = reshape(QPTIMETRANS(X, h, Nt), repmat(Nt, 1, Nc));
    xd  = reshape(QPTIMETRANS(D1*X, h, Nt), repmat(Nt, 1, Nc));
    xdd = reshape(QPTIMETRANS(D1^2*X, h, Nt), repmat(Nt, 1, Nc));
    
    x   = [x   x(:,1)  ; x(1,:)   x(1,1)];
    xd  = [xd  xd(:,1) ; xd(1,:)  xd(1,1)];
    xdd = [xdd xdd(:,1); xdd(1,:) xdd(1,1)];
    
    figure(4)
    clf()
    surf(x, xd, xdd, 'EdgeColor', 'r'); hold on
	plot3(U(1000:end), Ud(1000:end), Udd(1000:end), 'b-')
	grid on    
    xl = xlim; yl = ylim; 
    
    for ti=1:100:length(T)
        clf()
        % plot3(x, xd, xdd, '-', 'LineWidth', 1)
        surf(x, xd, xdd, 'EdgeColor', 'none', 'FaceAlpha', 1); hold on
        plot3(U(1:ti), Ud(1:ti), Udd(1:ti), '-', 'Color', [1 1 1]*0.6)
        plot3(U(ti), Ud(ti), Udd(ti), 'bo', 'MarkerFaceColor', 'b')
        grid on
        xlim(xl); ylim(yl);
        title(sprintf('Frame %d/%d', ti, length(T)))
        pause(0.01)
    end
end