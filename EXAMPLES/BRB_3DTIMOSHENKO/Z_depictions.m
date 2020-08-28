clc
clear all
addpath('../../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',14)

Nein = 8;
load(sprintf('./MATS/%dIN_MATS.mat',Nein), 'Beam1', 'Beam2', ...
    'SensorLocs', 'L1', 'L2', 'BoltLocs', 'IN1', 'IN2', 'wdt', ...
    'i1s', 'i2s', 'Nebb')

%% Bolt prestress
sig = (14.3e-3 + 2*tand(33)*25.4e-3)/6;
mus = BoltLocs;
fbolt = @(x, y) 1.0/(2*pi*(sig))*sum(exp(-0.5*((x-mus)/sig).^2-0.5*(y/sig).^2), 2);

%% Quadrature Points
No = 3;

[yi, wi] = LGWT(No, -wdt/2, wdt/2);
Ysi = repmat(kron(yi, ones(No, 1)), Nein, 1);

Les = diff(IN1.X);
Wys = IN1.WY(1:end-1);
Zs  = 0-IN1.Z(1:end-1);
[Q1, T1] = TM3D_ND2QP(Les, Wys, Zs, No);

Les = diff(IN2.X);
Wys = IN2.WY(1:end-1);
Zs  = 0-IN2.Z(1:end-1);
[Q2, T2] = TM3D_ND2QP(Les, Wys, Zs, No);

Q1p = zeros(Nein*No^2*3, (Nebb+Nein+1)*2*6);
Q2p = zeros(Nein*No^2*3, (Nebb+Nein+1)*2*6);

Q1p(:, i1s) = Q1;
Q2p(:, i1s) = Q2;

%% Plot
sis = 1:size(SensorLocs,1);

figure(1)
clf()
set(gcf, 'color', 'white')
set(gcf, 'Position', [1 1 1920 998])

DEPICTBEAM_TM3D(diff(Beam1.X), Beam1.WY, Beam1.WZ, ...
    [Beam1.X, Beam1.Y, Beam1.Z], zeros(size(L1(:,1))), 'b', ...
    0.1, 0);  % 0.1, 2
DEPICTBEAM_TM3D(diff(Beam2.X), Beam2.WY, Beam2.WZ, ...
	[Beam2.X, Beam2.Y, Beam2.Z], zeros(size(L2(:,1))), 'r', ...
	0.1, 0);

plot3(SensorLocs(sis,1), SensorLocs(sis,2), ...
    SensorLocs(sis,3), 'k.', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 12)
text(SensorLocs(sis,1), SensorLocs(sis,2), ...
    SensorLocs(sis,3), int2str(sis'))

axis equal
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
% ylim(25.4e-3*[-1 1])

axes('Position',[.1 .5 .4 .4])
Ng = 100;
[xgrid, ygrid] = meshgrid(linspace(IN1.X(1), IN1.X(end), Ng), linspace(-wdt/2, wdt/2, Ng));

surf(xgrid, ygrid, 1e-3*reshape(fbolt(xgrid(:), ygrid(:)), Ng, []), 'EdgeColor', 'None')
axis equal
grid off
box on
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])

xlabel('Spatial X', 'Rotation', 22)
ylabel('Spatial Y', 'Rotation', -35)
zlabel('Normal Traction')
title('Bolt Loading')

axes('Position',[.6 .1 .3 .3])
DEPICTBEAM_TM3D(diff(Beam1.X), Beam1.WY, Beam1.WZ, ...
    [Beam1.X, Beam1.Y, Beam1.Z], zeros(size(L1(:,1))), [185 185 255]/255, ...
    1, 0);  % 0.1, 2
DEPICTBEAM_TM3D(diff(Beam2.X), Beam2.WY, Beam2.WZ, ...
	[Beam2.X, Beam2.Y, Beam2.Z+25.4e-3*2/2], zeros(size(L2(:,1))), [255 185 185]/255, ...
	1, 0);

plot3(Q1(1:3:end, 1:6:end)*IN1.X, Ysi, ones(size(Ysi))*0, 'k.')
plot3(Q2(1:3:end, 1:6:end)*IN2.X, Ysi, ones(size(Ysi))*0+25.4e-3*2/2, 'k.')

axis equal
box on
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
title('Exploded View of Interface')

xlim([299.72-30 299.72+120+30]*1e-3)

% print(sprintf('./FIGS/%dMODEL_WSENS.eps', Nein), '-depsc')
% print(sprintf('./FIGS/%dMODEL_WSENS_o.svg', Nein), '-dsvg')
export_fig(sprintf('./FIGS/%dMODEL_WSENS.svg', Nein), '-dsvg')