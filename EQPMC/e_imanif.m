clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/export_fig/')

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if ~isOctave
  set(0,'defaultAxesTickLabelInterpreter', 'default');
  set(0,'defaultTextInterpreter','latex');
  set(0, 'DefaultLegendInterpreter', 'latex');
  set(0,'defaultAxesFontSize',13);
end

analyze = false;
plotout = true;
%% Parameters
M = eye(2);
K = [2 -1;-1 2];
C = 0.01*K;
b = 0.5;

MDL = MDOFGEN(M, K, C, eye(2));
MDL = MDL.SETNLFUN(1+3, [1 0], @(t, u, ud) deal(b*u.^3, 3*b*u.^2, ud*0));

[V, D] = eig(K, M);
[D, si] = sort(diag(D));
V = V(:,si);
%% EPMC
h = [0 1 2 3 4 5 6 7];
% h = 1:7;
Nt = 512;
Nhc = sum(h==0)+2*sum(h~=0);
Fl = zeros(Nhc*2,1);
Fl(sum(h==0)+2) = 1;

Astart = -2;
Aend = 4;
da = 0.1;
Copt = struct('Nmax', 2000, 'Display', 1, 'DynDscale', 1, 'dsmin', ...
              0.0001, 'angopt', 1e-1, 'dsmax', 0.5);
Copt.Dscale = [kron([1e-8; 1e-1*ones(Nhc-1,1)],[1;1]); 1.0; 0.01; 1.0];
% Copt.Dscale = [ones(Nhc*2,1)*1e-1; 1; 1e-2; 1.0];

BB = cell(2,1);
if analyze
    U0 = kron([zeros(h(1)==0,1); 1; 1; zeros(Nhc-2-sum(h==0),1)], V(:, 1));
    da = 0.001;
    BB{1} = CONTINUE(@(Uwxa) MDL.EPMCRESFUN(Uwxa, Fl, h, Nt, 1e-6), [U0; sqrt(D(1)); 0.01], Astart, Aend, da, Copt);
    
    U0 = kron([0; 1; 1; zeros(Nhc-3,1)], V(:, 2));
    BB{2} = CONTINUE(@(Uwxa) MDL.EPMCRESFUN(Uwxa, Fl, h, Nt, 1e-6), [U0; sqrt(D(2)); 0.01], Astart, Aend, da, Copt);
    
    save('./DATA/2dofdat.mat', 'BB');    
else
    load('./DATA/2dofdat.mat', 'BB');    
end

%% Invariant Manifold
D1 = HARMONICSTIFFNESS(0, 1, 0, 1, h);

for mi=1:2
    X = BB{mi}(1:2:end-3,:).*(10.^(BB{mi}(end,:)));
    Y = BB{mi}(2:2:end-3,:).*(10.^(BB{mi}(end,:)));
    Xdot = (D1*X).*BB{mi}(end-2,:);
    Ydot = (D1*Y).*BB{mi}(end-2,:);

    x = TIMESERIES_DERIV(Nt, h, X, 0);
    y = TIMESERIES_DERIV(Nt, h, Y, 0);
    xd = TIMESERIES_DERIV(Nt, h, Xdot, 0);
    yd = TIMESERIES_DERIV(Nt, h, Ydot, 0);

    ai = 95*mi;
    figure(mi)
    clf()
    surf(x([1:end 1], 1:ai), xd([1:end 1], 1:ai), y([1:end 1], 1:ai), 'EdgeColor', 'none'); hold on
    surf(x([1:end 1], 1:ai), xd([1:end 1], 1:ai), x([1:end 1], 1:ai)*V(1, mi)/V(2, mi), 'EdgeColor', 'none', 'FaceColor', [0 0 1]*0.6, 'FaceAlpha', 0.5)
    xlabel('Displacement $x_1(t)$')
    ylabel('Velocity $\dot{x_1}(t)$')
    zlabel('Displacement $x_2(t)$')
%     if mi==1
%         ylim([-1 1]*2)
%         caxis([-1 1]*20)
%         zlim([-1 1]*20)
%     else
%         xlim([-1 1]*6)
%         ylim([-1 1]*10);
%         zlim([-1 1]*1.2)
%         caxis([-1 1]*1.2)
%         set(gca, 'View', [20 25])
%     end
    colormap(jet)
%     print(sprintf('./FIGS/intro_2dofim%d.png',mi), '-dpng')
end

%% Detect "tongues"
[~, locs] = findpeaks(BB{1}(end,:)-min(BB{1}(end,:)));
[~, locis] = findpeaks(-BB{1}(end,:)+max(BB{1}(end,:)));

%% Frequency Amplitude Plot
% xlims = [10 800];
% xlims = [300 4e5];
xlims = [10 450];
ylims = [1.35 1.42];

figure(10)
clf()
semilogx(10.^BB{1}(end,:), BB{1}(end-2,:), '-', 'LineWidth', 2); hold on
semilogx(10.^BB{1}(end,locis), BB{1}(end-2,locis), 'o', 'LineWidth', 2); hold on

% semilogx(10.^BB{2}(end,:), BB{2}(end-2,:)/5, '-', 'LineWidth', 2)
xlabel('Modal Amplitude')
ylabel('Frequency (rad/s)')
xlim(10.^[Astart+0.1 Aend])
ylim([0.975 1.425])
plot(xlims([1 2 2 1 1]), ylims([1 1 2 2 1]), 'k-')

axes('Position', [.52 .2 .325 .325])
semilogx(10.^BB{1}(end,:), BB{1}(end-2,:), '-', 'LineWidth', 2); hold on
xlim(xlims)
ylim(ylims)
set(gca, 'XTick', xlims);
set(gca, 'YTick', ylims)

	      % print(sprintf('./FIGS/intro_2doffep.svg',mi), '-dsvg')

%% Plot frequency content
his = [2:2:h(end)];
figure(11)
clf()
fconts1 = sqrt(BB{1}([1 3:4:end-3],:).^2+[zeros(1,size(BB{1},2)); BB{1}([5:4:end-3],:).^2]);
fconts2 = sqrt(BB{1}([2 4:4:end-3],:).^2+[zeros(1,size(BB{1},2)); BB{1}([6:4:end-3],:).^2]);
loglog(10.^BB{1}(end,:), fconts1(his,:)./fconts1(2,:), '-', 'LineWidth', 2); hold on
loglog(10.^BB{1}(end,locis(1)), fconts1(his,locis(1))./fconts1(2,locis(1)), 'o', 'LineWidth', 2); hold on
% loglog(10.^BB{1}(end,:), fconts2(his,:)./fconts2(2,:), '.-', 'LineWidth', 2);
legend(arrayfun(@(c) sprintf('H%d', c), h(his), 'UniformOutput', false))

ylim([1e-3 3e1])
xlabel('Modal Amplitude')
ylabel('Relative Harmonic Content in x1')

%% Plot Parseval's Amplitude
mi = 1;

figure(12)
clf()
fpars = sqrt(([0 0.5*ones(1, Nhc-1)]*BB{mi}(1:2:end-3,:).*(10.^BB{mi}(end,:))).^2);
fh1 = sqrt(([0 0.5*ones(1, 2) zeros(1, Nhc-3)]*BB{mi}(1:2:end-3,:).*(10.^BB{mi}(end,:))).^2);
plot3(fh1, fpars, BB{mi}(end-2,:), '-', 'LineWidth', 2)
grid on
set(gca, 'Xscale', 'log', 'Yscale', 'log')
xlabel('$|H1|$')
ylabel('$|RMS|$')
zlabel('Natural Frequency')

%% Plot Harmonic Mode Shape Contributions
MS1 = reshape(BB{1}(3:end-3, :), 2, Nhc-1, size(BB{1},2));
MS1 = MS1(:, 1:2:end, :)-1j*MS1(:, 2:2:end, :);
% MSn1 = reshape(sqrt(sum((M*reshape(MS1, 2, [])).*reshape(MS1, 2, []), 1))', length(h)-1, length(EN1));

MACs_11 = squeeze(sum((V(:, 1)'*M)'.*MS1, 1));
MACs_21 = squeeze(sum((V(:, 2)'*M)'.*MS1, 1));

figure(13)
clf()
semilogx(10.^BB{1}(end,:), abs(MACs_21([1 3 5],:)))
grid on
xlabel('Modal Amplitude')
ylabel('$MAC_{21}$')

%% Correlate Interactions by Energy
% KE1 = ((10.^BB{1}(end,:)).*BB{1}(end-2,:)).^2;
% KE2 = ((10.^BB{2}(end,:)).*BB{2}(end-2,:)).^2;

Uh1 = BB{1}(3:end-3,:).*(10.^BB{1}(end,:));
Uh2 = BB{2}(3:end-3,:).*(10.^BB{2}(end,:));

Vh1 = Uh1.*BB{1}(end-2,:);
Vh2 = Uh2.*BB{2}(end-2,:);

KE1 = 1/2*0.5*sum((kron(eye(Nhc-1),M)*Vh1).*Vh1,1);
KE2 = 1/2*0.5*sum((kron(eye(Nhc-1),M)*Vh2).*Vh2,1);

PE1 = 1/2*0.5*sum((kron(eye(Nhc-1),K)*Uh1).*Uh1,1) + 1/4*3.0/8.0*b*sum(Uh1(1:2:end,:).^4,1);
PE2 = 1/2*0.5*sum((kron(eye(Nhc-1),K)*Uh2).*Uh2,1) + 1/4*3.0/8.0*b*sum(Uh2(1:2:end,:).^4,1);

EN1 = KE1+PE1;
EN2 = KE2+PE2;

Npc = 100;
comgrid = EN1(find(EN1>=EN2(1), 1):find(EN1<=EN2(end), 1, 'last'));
bb1 = interp1(EN1, BB{1}', comgrid)';
bb2 = interp1(EN2, BB{2}', comgrid)';

[~, locis] = findpeaks(-comgrid+max(comgrid));

figure(14)
clf()
semilogx(comgrid, bb1(end-2,:), '-'); hold on
semilogx(comgrid, bb2(end-2,:)/3, '-')

plot(comgrid(locis(1))*[1 1], ylim, 'k--')
plot(comgrid(locis(1)), bb1(end-2,locis(1)), 'ko')
plot(comgrid(locis(1)), bb2(end-2,locis(1))/3, 'k*')
ylim([0.975 1.425])
xlabel('Energy (J)')
ylabel('Natural Frequency')
% semilogx(EN2, BB{2}(end-2,:)/3, '-'); hold on

%% Correlate Mode Shapes on Common Grid Points
MS1 = reshape(bb1(3:end-3,:), 2, Nhc-1, length(comgrid));
MS1 = MS1(:, 1:2:end, :)-1j*MS1(:, 2:2:end, :);
MSn1 = reshape(sqrt(sum((M*reshape(MS1, 2, [])).*reshape(MS1, 2, []), 1))', length(h)-1, length(comgrid));
% MSn1 = 1;

MS2 = reshape(bb2(3:end-3,:), 2, Nhc-1, length(comgrid));
MS2 = MS2(:, 1:2:end, :)-1j*MS2(:, 2:2:end, :);
MSn2 = reshape(sqrt(sum((M*reshape(MS2, 2, [])).*reshape(MS2, 2, []), 1))', length(h)-1, length(comgrid));
% MSn2 = 1;

MACs_21 = zeros(length(h)-(h(1)==0), length(comgrid));
MACs_12 = zeros(length(h)-(h(1)==0), length(comgrid));

it = 1;
MACs_21 = squeeze(sum(permute(squeeze(MS2(:, 1, [it:end 1:it-1]))'*M, [2 3 1]).*MS1, 1))./MSn1;
MACs_12 = squeeze(sum(permute(squeeze(MS1(:, 1, [it:end 1:it-1]))'*M, [2 3 1]).*MS2, 1))./MSn2;

figure(15)
clf()
semilogx(comgrid, abs(MACs_21(1:2:end,:)), '-'); hold on
semilogx(comgrid(locis(1)), abs(MACs_21(1:2:end,locis(1))), 'ko')
% semilogx(comgrid, abs(MACs_12(1:2:end,:)), '-')
xlabel('Energy (J)')
ylabel('$MAC_{21}$')
title(sprintf('Max = %f\n', max(abs(MACs_21(3,:)))))

% !!!

%% Plot mode shapes in time
ms1 = TIMESERIES_DERIV(Nt, h, reshape(permute(reshape(BB{1}(1:end-3,:), 2, size(BB{1},2), Nhc), [3 1 2]), Nhc, []), 0);
ms1 = reshape(ms1, Nt, 2, []);
ms2 = TIMESERIES_DERIV(Nt, h, reshape(permute(reshape(BB{2}(1:end-3,:), 2, size(BB{2},2), Nhc), [3 1 2]), Nhc, []), 0);
ms2 = reshape(ms2, Nt, 2, []);

figure(200)
clf()
plot(ms1(:, 1:2, locis(1)));
