clc
clear all
addpath('../../ROUTINES/FEM')
addpath('../../ROUTINES/FEM/BEAMS')

%% Set Geometric Properties
Lbb = 11.8*25.4e-3;  % Body length
Lin = 120e-3;  % Interface Length
wdt = 25.4e-3;

BoltLocs = Lbb+30e-3*(1:3);  % Location of bolts

%% Set Material Properties
nu = 0.3;  % Poisson's ratio
kappa = 5/6;  % Traditional
kappa = 10*(1+nu)/(12+11*nu);  % Cowper (1966) for rectangular sections
pars = struct('E', 2e11, 'rho', 7800, 'G', 2e11/(2*(1+nu)), ...
	      'A', wdt*wdt, 'Iy', wdt^4/12, 'Iz', wdt^4/12, ...
	      'k1', kappa, 'k2', kappa, 'k3', kappa);

parsint = pars;  % Half beam portion
parsint.A = pars.A/2;
parsint.Iy = pars.Iy/8;
parsint.Iz = pars.Iz/2;

%% Set Finite Element Discretization
Nein = 8;  % Number of elements in the interface
Nebb = fix(Nein*2.5);  % Number of elements in the rest of the body

% Coordinates of nodes (neutral axis)
BM1 = struct('X', linspace(0, Lbb, Nebb+1)', ...
	     'Y', zeros(Nebb+1,1), ...
	     'Z', zeros(Nebb+1,1), ...
             'WY', wdt*ones(Nebb+1,1), ...
             'WZ', wdt*ones(Nebb+1,1));
IN1 = struct('X', linspace(Lbb, Lbb+Lin, Nein+1)', ...
	     'Y', zeros(Nein+1,1), ...
	     'Z', ones(Nein+1,1)*(-wdt/4), ...
             'WY', wdt*ones(Nein+1,1), ...
             'WZ', wdt/2*ones(Nein+1,1));
Beam1 = struct('X', [BM1.X; IN1.X], ...
	       'Y', [BM1.Y; IN1.Y], ...
	       'Z', [BM1.Z; IN1.Z], ...
               'WY', [BM1.WY; IN1.WY], ...
               'WZ', [BM1.WZ(1:end-1); 0; IN1.WZ]);

IN2 = struct('X', linspace(Lbb, Lbb+Lin, Nein+1)', ...
	     'Y', zeros(Nein+1,1), ...
	     'Z', ones(Nein+1,1)*wdt/4, ...
             'WY', wdt*ones(Nein+1,1), ...
             'WZ', wdt/2*ones(Nein+1,1));
BM2 = struct('X', linspace(Lbb+Lin, Lbb+Lin+Lbb, Nebb+1)', ...
	     'Y', zeros(Nebb+1,1), ...
	     'Z', zeros(Nebb+1,1), ...
             'WY', wdt*ones(Nebb+1,1), ...
             'WZ', wdt*ones(Nebb+1,1));
Beam2 = struct('X', [IN2.X; BM2.X], ...
	       'Y', [IN2.Y; BM2.Y], ...
	       'Z', [IN2.Z; BM2.Z], ...
               'WY', [IN2.WY; BM2.WY], ...
	       'WZ', [IN2.WZ(1:end-1); 0; BM2.WZ]);

%% Construct Matrices

%% BOLT IDEALIZATION
%% BOLT ELEMENTS
rad = 8.42e-3/2;  % Bolt shank radius
kappa = 6*(1+nu)/(7+6*nu);

parsbolt = struct('E', 2e11, 'rho', 7800, 'G', 0*2e11/(2*(1+nu)), ...
		  'A', pi*rad^2, 'Iy', pi*rad^4/4, 'Iz', pi*rad^4/4, ...
		  'k1', kappa, 'k2', kappa, 'k3', kappa);
[Mbolt, Kbolt] = TM3D_ELMATS(25.4e-3, parsbolt);

mhead = 2.864e-2 - sum(diag(Mbolt(1:6:end,1:6:end)))*6/4;  % Discrete mass at head and nut
Mbolt([1 2 3 7 8 9], [1 2 3 7 8 9]) = mhead/2*eye(6);
Mbolt([4 5 6 10 11 12], [4 5 6 10 11 12]) = mhead*25.4e-3^2*eye(6);

km = 1.2041e9;  % Member stiffness due to bolt prestress
Kbolt([1 7], [1 7]) = Kbolt([1 7], [1 7]) + km*[1 -1;-1 1];

% Kbolt([2 8], [2 8]) = Kbolt([1 7], [1 7]) + km*[1 -1;-1 1];  % Not possible

%% Bolt CS (X, Y, Z) => Global CS (Z, X, Y)
TFM = [0 0 -1;0 1 0;1 0 0];
TFM = kron(eye(4), TFM);

Mbolt = TFM'*Mbolt*TFM;
Kbolt = TFM'*Kbolt*TFM;

% Bolt Nodes (on the interface
Bnodes = fix(Nein*2.5);
if Nein>=4
  Bnodes = Nein/4*[1 2 3]+1;
end

%% Bolt prestress shape based on gaussian pressure distributions
%% set up quadrature
No = 20;
[yi, wi] = LGWT(No, -wdt/2, wdt/2);

Les = diff(IN1.X);
Wys = IN1.WY(1:end-1);
Zs  = -wdt/2-IN1.Z(1:end-1);
[Q1, T1] = TM3D_ND2QP(Les, Wys, Zs, No);  % Bottom Face

Les = diff(IN2.X);
Wys = IN2.WY(1:end-1);
Zs  = wdt/2-IN2.Z(1:end-1);
[Q2, T2] = TM3D_ND2QP(Les, Wys, Zs, No);  % Top Face

sig = (14.3e-3 + 2*tand(33)*25.4e-3)/6;
mus = BoltLocs;
fbolt = @(x, y) 1.0/(2*pi*(sig))*sum(exp(-0.5*((x-mus)/sig).^2-0.5*(y/sig).^2), 2);

Ysi = repmat(kron(yi, ones(No, 1)), Nein, 1);

Xsi = Q1(1:3:end, 1:6:end)*IN1.X;
Fb_1 = T1(:, 3:3:end)*fbolt(Xsi, Ysi);
fscal = sum(Fb_1(3:6:end));
Fb_1 = 3*Fb_1/fscal;  % Normalize sum to apply 3N 

Xsi = Q2(1:3:end, 1:6:end)*IN2.X;
Fb_2 = T2(:, 3:3:end)*fbolt(Xsi, Ysi);
Fb_2 = 3*Fb_2/sum(Fb_2(3:6:end));  % Normalize sum to apply 3N 

% Ng = 100;
% [xgrid, ygrid] = meshgrid(linspace(IN1.X(1), IN1.X(end), Ng), linspace(-wdt/2, wdt/2, Ng));

% figure(1)
% clf()
% surf(xgrid, ygrid, fscal*reshape(fbolt(xgrid(:), ygrid(:)), Ng, []), 'EdgeColor', 'None'); hold on

% # plot3(IN1.X, IN1.Y, Fb_1(3:6:end), 'o');
% for ni=1:length(IN1.X)
%   plot3(IN1.X(ni)*[1 1], IN1.Y(ni)*[1 1], Fb_1((ni-1)*6+3)*[0 1], 'ko-')
% end
% xlabel('X Coordinate')
% ylabel('Y Coordinate')
% zlabel('Forcing')

% figure(2)
% clf()
% plot(Xsi, Ysi, '*'); hold on
% % text(Xsi, Ysi, int2str((1:Nein*No^2)'))
% plot(IN1.X, IN1.X*0, 'ko-')
% grid on

% xlim([IN1.X(1) IN1.X(end)])
% ylim(wdt*[-1 1])

%% BEAM 1
% Body 1
[M1, K1] = CONSTRUCT(BM1, pars);
% Interface 1
[Mi1, Ki1] = CONSTRUCT(IN1, parsint);
% Tie the last node of BM1 to 1st node of IN1 by placing IN1 on BM1's section
L1 = eye((Nebb+1+Nein+1)*6);
bnode = Nebb+1;  % Last node of body
bis   = (bnode-1)*6+(1:6);
inode = Nebb+1+1;  % First node in interface
eis   = (inode-1)*6+(1:6);

dz = BM1.Z(1)-IN1.Z(end);
dy = BM1.Y(1)-IN1.Y(end);
L1(bis, eis) = [eye(3), [0 dz -dy; -dz 0 0; dy 0 0];
		zeros(3), eye(3)];  % Represent body node as a function of interface node
L1(:, bis) = [];   % Only retain node on the interface

Mbm1 = L1'*blkdiag(M1, Mi1)*L1;
Kbm1 = L1'*blkdiag(K1, Ki1)*L1;

%% BEAM 2
% Interface 2
[Mi2, Ki2] = CONSTRUCT(IN2, parsint);
% Body 2
[M2, K2] = CONSTRUCT(BM2, pars);
% Tie the last node of BM1 to 1st node of IN1 by placing IN1 on BM1's section
L2 = eye((Nein+1+Nebb+1)*6);
inode = Nein+1;  % Last node of interface
eis   = (inode-1)*6+(1:6);
bnode = Nein+1+1;  % First node of body
bis   = (bnode-1)*6+(1:6);

dz = BM2.Z(1)-IN2.Z(end);
dy = BM2.Y(1)-IN2.Y(end);
L2(bis, eis) = [eye(3), [0 dz -dy; -dz 0 0; dy 0 0];
		zeros(3), eye(3)];   % Represent body node as a function of interface node
L2(:, bis) = [];   % Retain node on the interface

Mbm2 = L2'*blkdiag(Mi2, M2)*L2;
Kbm2 = L2'*blkdiag(Ki2, K2)*L2;

%% Recovery Matrix for sensor/actuator locations
sens1 = [0; 100; 200]*1e-3;
[ndi, ~] = find((Beam1.X(1:end-1)-sens1').*(Beam1.X(2:end)-sens1')<=0);
SFs = [(Beam1.X(ndi+1)-sens1)./(Beam1.X(ndi+1)-Beam1.X(ndi)), ...
       (Beam1.X(ndi)-sens1)./(Beam1.X(ndi)-Beam1.X(ndi+1))];

R1 = zeros(length(sens1), length(Beam1.X));
R1(:, [ndi ndi+1]) = [diag(SFs(:,1)) diag(SFs(:,2))];

don = [eye(3) zeros(3)];  % Displacements at node
zpd = [eye(3) [0 -wdt/2 0;wdt/2 0 0;0 0 0]];  % Displacement at z+
zmd = [eye(3) [0 wdt/2 0;-wdt/2 0 0;0 0 0]];  % Displacement at z-
ypd = [eye(3) [0 0 wdt/2;0 0 0;-wdt/2 0 0]];  % Displacement at y+
ymd = [eye(3) [0 0 -wdt/2;0 0 0;wdt/2 0 0]];  % Displacement at y+

% R1 = kron(R1, [eye(3) zeros(3)])*L1*[eye((Nebb+Nein+1)*6) zeros((Nebb+Nein+1)*6)];  % Only displacements at Nodes

R1 = [kron(R1(1,:), don); kron(R1(2:end,:), [zpd; ypd; zmd; ymd])]*L1*[eye((Nebb+Nein+1)*6) zeros((Nebb+Nein+1)*6)];

sens2 = Beam2.X(end) - [200; 100; 0]*1e-3;
[ndi, ~] = find((Beam2.X(1:end-1)-sens2').*(Beam2.X(2:end)-sens2')<=0);

SFs = [(Beam2.X(ndi+1)-sens2)./(Beam2.X(ndi+1)-Beam2.X(ndi)), ...
       (Beam2.X(ndi)-sens2)./(Beam2.X(ndi)-Beam2.X(ndi+1))];

R2 = zeros(length(sens2), length(Beam2.X));
R2(:, [ndi ndi+1]) = [diag(SFs(:,1)) diag(SFs(:,2))];

don = [eye(3) zeros(3)];  % Displacements at node
zpd = [eye(3) [0 -wdt/2 0;wdt/2 0 0;0 0 0]];  % Displacement at z+
zmd = [eye(3) [0 wdt/2 0;-wdt/2 0 0;0 0 0]];  % Displacement at z-
ypd = [eye(3) [0 0 wdt/2;0 0 0;-wdt/2 0 0]];  % Displacement at y+
ymd = [eye(3) [0 0 -wdt/2;0 0 0;wdt/2 0 0]];  % Displacement at y+
% R2 = kron(R2, [eye(3) zeros(3)])*L2*[zeros((Nebb+Nein+1)*6) eye((Nebb+Nein+1)*6)];

R2 = [kron(R2(1:end-1,:), [zpd; ypd; zmd; ymd]); kron(R2(end,:), don)]*L2*[zeros((Nebb+Nein+1)*6) eye((Nebb+Nein+1)*6)];

RECOV = [R1; R2];

SensorLocs = [sens1(1) 0 0;   % 1st node from left
	    %
	    sens1(2) 0 wdt/2;    % z+ 100mm from left
	    sens1(2) wdt/2 0;    % y+ 100mm from left
	    sens1(2) 0 -wdt/2;   % z- 100mm from left
	    sens1(2) -wdt/2 0;   % y- 100mm from left
	    %
	    sens1(3) 0 wdt/2;    % z+ 200mm from left
	    sens1(3) wdt/2 0;    % y+ 200mm from left
	    sens1(3) 0 -wdt/2;   % z- 200mm from left
	    sens1(3) -wdt/2 0;   % y- 200mm from left
	    %
	    sens2(1) 0 wdt/2;    % z+ 200mm from right
	    sens2(1) wdt/2 0;    % y+ 200mm from right
	    sens2(1) 0 -wdt/2;   % z- 200mm from right
	    sens2(1) -wdt/2 0;   % y- 200mm from right
	    %
	    sens2(2) 0 wdt/2;    % z+ 100mm from right
	    sens2(2) wdt/2 0;    % y+ 100mm from right
	    sens2(2) 0 -wdt/2;   % z- 100mm from right
	    sens2(2) -wdt/2 0;   % y- 100mm from right
	    %
	    sens2(3) 0 0];   % 1st node from right

%% Full System Assembly
M = blkdiag(Mbm1, Mbm2);
K = blkdiag(Kbm1, Kbm2);

%% Bolts
Bnodes1 = Nebb+Bnodes;
Bnodes2 = Nebb+Nein+1 + Bnodes;

for i = 1:3
  K([(Bnodes1(i)-1)*6+(1:6) (Bnodes2(i)-1)*6+(1:6)], ...
    [(Bnodes1(i)-1)*6+(1:6) (Bnodes2(i)-1)*6+(1:6)]) = ...
  K([(Bnodes1(i)-1)*6+(1:6) (Bnodes2(i)-1)*6+(1:6)], ...
    [(Bnodes1(i)-1)*6+(1:6) (Bnodes2(i)-1)*6+(1:6)]) + Kbolt;

  M([(Bnodes1(i)-1)*6+(1:6) (Bnodes2(i)-1)*6+(1:6)], ...
    [(Bnodes1(i)-1)*6+(1:6) (Bnodes2(i)-1)*6+(1:6)]) = ...
  M([(Bnodes1(i)-1)*6+(1:6) (Bnodes2(i)-1)*6+(1:6)], ...
    [(Bnodes1(i)-1)*6+(1:6) (Bnodes2(i)-1)*6+(1:6)]) + Mbolt;
end

int1nds = Nebb+(1:(Nein+1));
int2nds = (Nebb+Nein+1)+(1:(Nein+1));
remnds  = setdiff(1:(Nebb+Nein+1)*2, [int1nds, int2nds]);

i1s = reshape((int1nds-1)*6+(1:6)', [], 1);
i2s = reshape((int2nds-1)*6+(1:6)', [], 1);

Fbolt = zeros(size(K,1),1);
Fbolt([i1s; i2s]) = [Fb_1; -Fb_2];

%% Analytical Rigid Body Modes (in the stuck state)
Xcoords = [BM1.X(1:end-1); IN1.X; IN2.X; BM2.X(2:end)];
Ycoords = [BM1.Y(1:end-1); IN1.Y; IN2.Y; BM2.Y(2:end)];
Zcoords = [BM1.Z(1:end-1); IN1.Z; IN2.Z; BM2.Z(2:end)];

Vrbms = zeros(size(K,1), 6);
Vrbms(1:6:end,1) = 1;
Vrbms(2:6:end,2) = 1;
Vrbms(3:6:end,3) = 1;

Vrbms(4:6:end,4) = 1;
Vrbms(2:6:end,4) = -Zcoords;
Vrbms(3:6:end,4) = Ycoords;

Vrbms(5:6:end,5) = 1;
Vrbms(1:6:end,5) = Zcoords;
Vrbms(3:6:end,5) = -Xcoords;

Vrbms(6:6:end,6) = 1;
Vrbms(1:6:end,6) = -Ycoords;
Vrbms(2:6:end,6) = Xcoords;

%% Save the Necessities
L1 = sparse(L1*[eye((Nebb+Nein+1)*6) zeros((Nebb+Nein+1)*6)]);
L2 = sparse(L2*[zeros((Nebb+Nein+1)*6) eye((Nebb+Nein+1)*6)]);
M = sparse(M); K = sparse(K); Fbolt = sparse(Fbolt); Vrbms = sparse(Vrbms);
R1 = sparse(R1); R2 = sparse(R2);
RECOV = sparse(RECOV);

save(sprintf('./MATS/%dIN_MATS.mat',Nein), 'M', 'K', 'Fbolt', 'L1', 'L2', ...
     'SensorLocs', 'RECOV', 'R1', 'R2', 'BM1', 'IN1', 'Beam1', ...
     'BM2', 'IN2', 'Beam2', 'pars', 'parsint', 'parsbolt', 'Nebb', ...
     'Nein', 'wdt', 'nu', 'int1nds', 'int2nds', 'remnds', ...
     'i1s', 'i2s', 'Xcoords', 'Ycoords', 'Zcoords', 'Vrbms', '-v7');

% Check Saved Data
clear all
Nein = 8;
load(sprintf('./DATA/%dIN_MATS.mat',Nein))

%% SET UP QUADRATURE
No = 2;

Les = diff(IN1.X);
Wys = IN1.WY(1:end-1);
Zs  = 0-IN1.Z(1:end-1);
[Q1, T1] = TM3D_ND2QP(Les, Wys, Zs, No);

Les = diff(IN2.X);
Wys = IN2.WY(1:end-1);
Zs  = 0-IN2.Z(1:end-1);
[Q2, T2] = TM3D_ND2QP(Les, Wys, Zs, No);

% Assembly for relative coordinates at quadrature locations
int1nds = Nebb+(1:(Nein+1));
int2nds = (Nebb+Nein+1)+(1:(Nein+1));
remnds  = setdiff(1:(Nebb+Nein+1)*2, [int1nds, int2nds]);

i1s = reshape((int1nds-1)*6+(1:6)', [], 1);
i2s = reshape((int2nds-1)*6+(1:6)', [], 1);

%% IN1 is facing downwards (z-) and IN2 is facing upwards (z+) so
%% positive interference is detected by IN2_disp-IN1_disp and this is
%% the sense in which relative displacements will be interpreted.

Qrel = zeros(Nein*No^2*3, (Nebb+Nein+1)*2*6);
Trel = zeros((Nebb+Nein+1)*2*6, Nein*No^2*3);

Qrel(:, i1s) = -Q1;
Qrel(:, i2s) = Q2;

Trel(i1s, :) = -T1;
Trel(i2s, :) = T2;

%% Contact Model (Hertz-Mindlin)
Prestress = 12e3;
Aint = sum(sum(T1(1:6:end, :)));

Pint = Prestress/Aint;
sint = 15e-6;
chi  = 2.0;
ktkn = chi*(1-nu)/(2-nu);
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn   = kt/ktkn;
% kn = 1e10;

tstiff = zeros(Nein*No^2*3,1);
tstiff(1:3:end) = kt;
tstiff(2:3:end) = kt;
tstiff(3:3:end) = kn;

Ks = Trel*diag(tstiff)*Qrel;

[V, W] = eig(full(K), full(M));
[W, si] = sort(sqrt(diag(W))/2/pi);
V = V(:,si);
V = V./sqrt(diag(V'*M*V))';

sensdat = RECOV*V;

[V0, W0] = eig(full(K+Ks), full(M));
[W0, si] = sort(sqrt(diag(W0))/2/pi);
V0 = V0(:,si);
V0 = V0./sqrt(diag(V0'*M*V0))';

sensdat = RECOV*V0;

%% Plot Mode Shapes

mi = 7;
sc = 1e-1;

figure(1)
clf()
DEPICTBEAM_TM3D(diff(Beam1.X), Beam1.WY, Beam1.WZ, ...
		[Beam1.X, Beam1.Y, Beam1.Z], L1*V0(:, mi)*sc, 'b', ...
		0.1, 2);  % 0.1, 2
DEPICTBEAM_TM3D(diff(Beam2.X), Beam2.WY, Beam2.WZ, ...
		[Beam2.X, Beam2.Y, Beam2.Z], L2*V0(:, mi)*sc, 'r', ...
		0.1, 0);
axis equal
grid on

% plot3(Beam1.X+sc*L1(1:6:end,:)*V0(:, mi), Beam1.Y+sc*L1(2:6:end,:)*V0(:, mi), ...
%       Beam1.Z+sc*L1(3:6:end,:)*V0(:, mi), 'bo-'); hold on
% plot3(Beam2.X+sc*L2(1:6:end,:)*V0(:, mi), Beam2.Y+sc*L2(2:6:end,:)*V0(:, mi), ...
%       Beam2.Z+sc*L2(3:6:end,:)*V0(:, mi), 'ro-'); hold on
plot3(SensorLocs(:,1)+sc*sensdat(1:3:end, mi), SensorLocs(:,2)+sc*sensdat(2:3:end, mi), ...
      SensorLocs(:,3)+sc*sensdat(3:3:end, mi), 'k.', 'MarkerFaceColor', 'k')

text(SensorLocs(:,1)+sc*sensdat(1:3:end, mi), SensorLocs(:,2)+sc*sensdat(2:3:end, mi), ...
     SensorLocs(:,3)+sc*sensdat(3:3:end, mi), int2str((1:size(SensorLocs,1))'))
     
title(sprintf('Mode %d: %.2f Hz', mi, W0(mi)))

xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')

disp(W0(1:10))

disp(diag(Vrbms'*(K+Ks)*Vrbms))

%% Linear Prestress Test
L = null(full(Vrbms'*M));
Ups = L*((L'*(K+Ks)*L)\(L'*Fbolt*12e3));
sc = 1.0

figure(1)
clf()
DEPICTBEAM_TM3D(diff(Beam1.X), Beam1.WY, Beam1.WZ, ...
		[Beam1.X, Beam1.Y, Beam1.Z], L1*Ups*sc, 'b', ...
		0.1, 2);  % 0.1, 2
DEPICTBEAM_TM3D(diff(Beam2.X), Beam2.WY, Beam2.WZ, ...
		[Beam2.X, Beam2.Y, Beam2.Z], L2*Ups*sc, 'r', ...
		0.1, 0);
% axis equal
grid on
ylim(25.4e-3*[-1 1])
xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')

figure(2)
clf()

plot(Beam1.X, L1(3:6:end,:)*Ups, 'bo'); hold on
plot(Beam2.X, L2(3:6:end,:)*Ups, 'ro')

plot(Q1(1:3:end,1:6:end)*IN1.X, Q1(3:3:end,:)*Ups(i1s), 'b-')
plot(Q2(1:3:end,1:6:end)*IN2.X, Q2(3:3:end,:)*Ups(i2s), 'r-')

grid on
xlabel('X Coordinate (m)')
ylabel('Z Displacement (m)')
