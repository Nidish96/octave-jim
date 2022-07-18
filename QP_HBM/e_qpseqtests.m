% testing different versions of hysteretic force evaluation for elastic dry
% friction
clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

figsaves = true;
%% Frequency configuration
Nc = 2;
ws = [1.0 sqrt(2)];
% ws = [1.0 pi];
% ws = [pi 1.0];

%%
Nt = 8;

t = linspace(0, 2*pi, Nt+1);  t(end) = [];
taus = cell(Nc, 1);
[taus{:}] = ndgrid(t);
t = t';

%% ranking based on "d" - No physical need for this
ds = (taus{1}*ws(1)+taus{2}*ws(2))/vecnorm(ws);
[~, si] = sort(ds(:));
% (i-1)*Nt+(j-1)

sxi = floor(si/Nt)+1;
sxi(sxi>Nt) = Nt;
syi = si - (sxi-1)*Nt;
sxi(syi==0) = sxi(syi==0)-1;
syi = si - (sxi-1)*Nt;

%% "Previous" points in time
delt = t(2)/vecnorm(ws);

figure(2)
clf()
for i=0:2
    for j=0:2
%         surf(taus{1}+(i-1)*t(end), taus{2}+(j-1)*t(end), ds*0, ...
%             'EdgeColor', [1 1 1]*0.6*(i~=1 | j~=1), 'FaceColor', 'none')
        surf(taus{1}+(i-1)*2*pi, taus{2}+(j-1)*2*pi, reshape(-(1:Nt^2), Nt,Nt), ...
            'EdgeColor', [1 1 1]*0.6*(i~=1 | j~=1)); hold on
    end
end
surf(taus{1}, taus{2}, ds*0, ...
    'EdgeColor', 'k', 'FaceColor', 'none'); hold on
set(gca, 'View', [0 90])
colormap(jet((Nt+1)^2))

angs = linspace(pi, 3*pi/2, 45);
for i=1:Nt
    for j=1:Nt
         plot(taus{1}(i,j) + t(2)*cos(angs), taus{2}(i,j) + t(2)*sin(angs), 'm-', 'LineWidth', 2)
    end
end

quiver(taus{1}-ws(1)*delt, taus{2}-ws(2)*delt, ones(Nt)*ws(1)*delt, ones(Nt)*ws(2)*delt, 'm-', 'LineWidth', 2)
plot(taus{1}, taus{2}, 'ko', 'MarkerFaceColor', 'w', 'LineWidth', 2); hold on
plot(taus{1}-ws(1)*delt, taus{2}-ws(2)*delt, 'ro', 'MarkerFaceColor', 'w', 'LineWidth', 2)
axis equal

xlim([-t(2)*2.25 2*pi+t(2)*1.5])
ylim([-t(2)*2.25 2*pi+t(2)*1.5])

set(gca, 'XTick', -2*pi:t(2):4*pi)
set(gca, 'YTick', -2*pi:t(2):4*pi)

%%
figure(1)
clf()
set(gcf, 'Color', 'white')
plot([-t(2) 2*pi], t.*[1 1], 'Color', [1 1 1]*.6, 'LineWidth', 0.25); hold on
plot(t.*[1 1], [-t(2) 2*pi], 'Color', [1 1 1]*.6, 'LineWidth', 0.25)
mesh(taus{1}, taus{2}, taus{1}*0, 'EdgeColor', 'k', 'LineWidth', 1); hold on
plot(taus{1}, taus{2}, 'ko', 'MarkerFaceColor', 'w'); hold on
axis equal

Np = 8;
am1 = 20;
am2 = 3.5;
aa = gobjects(2);
cols = repmat([1 0 0], 2*Np+1, 1);
kk = 1;

xsh = pi+t(2)-(2*pi*3/Np+(ws(1)/ws(2))*(pi+t(2)));
for i=-Np:Np
    quiver(2*pi*i/Np+xsh, -t(2), ws(1)/vecnorm(ws)*am1, ws(2)/vecnorm(ws)*am1, '-', 'LineWidth', 1, 'Color', cols(kk,:));
    quiver(2*pi*i/Np+xsh, -t(2), ws(1)/vecnorm(ws)*am2, ws(2)/vecnorm(ws)*am2, '-', 'LineWidth', 1.5, 'Color', cols(kk,:));
    kk = kk+1;
end

tt = text(taus{1}(:)-0.4, taus{2}(:)-0.2, num2str((1:Nt^Nc)'), 'FontWeight', 'bold');

tics = (0:.25:1);
ticls = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
set(gca, 'XTick', tics*2*pi)
set(gca, 'XTickLabel', ticls)
set(gca, 'YTick', tics*2*pi)
set(gca, 'YTickLabel', ticls)

xlim([-t(2) 2*pi])
ylim([-t(2) 2*pi])
xlabel('Scaled Time $\tau_1$')
ylabel('Scaled Time $\tau_2$')

pti = 37;
plot(taus{1}(pti)-t(2)*cos([t;2*pi]/4), taus{2}(pti)-t(2)*sin([t;2*pi]/4), 'b-', 'LineWidth', 1.5)
mm = plot(taus{1}(pti)-ws(1)*t(2)/vecnorm(ws), taus{2}(pti)-ws(2)*t(2)/vecnorm(ws), 'kh', 'MarkerFacecolor', 'b');
rr = plot(taus{1}(pti), taus{2}(pti), 'kh', 'MarkerFaceColor', 'r') ;

mesh(taus{1}(pti)+[-1.5 -1.5;.5 .5]*t(2), taus{2}(pti)+[-1.5 .5;-1.5 .5]*t(2), zeros(2), 'FaceColor', 'none', 'LineStyle', '-.', 'EdgeColor', [0 1 0]*0.6, 'LineWidth', 2)

if figsaves
    export_fig('./FIGS/E_qpgrid.png', '-dpng')
end
% %%
xl = taus{1}(pti)+[-1.5 .5]*t(2);
yl = taus{2}(pti)+[-1.5 .5]*t(2);
xlim(xl); ylim(yl);
for i=1:length(tt)
    if ~any(i==[pti pti-1 pti-Nt pti-Nt-1])
        tt(i).Visible = 'off';
    else
        tt(i).Position = tt(i).Position+[0.3 0.25 0];
        tt(i).FontSize = 15;
    end
end
quiver(taus{1}(pti)-ws(1)*t(2)/vecnorm(ws), taus{2}(pti)-ws(2)*t(2)/vecnorm(ws), ws(1)*t(2)/vecnorm(ws), ws(2)*t(2)/vecnorm(ws), '-', 'LineWidth', 2, 'Color', [1 1 0]*0.5)
rr.MarkerSize = 15;
mm.MarkerSize = 15;

if figsaves
    export_fig('./FIGS/E_qpgrid_zoom.png', '-dpng')
end

%%
figure(10)
clf()
set(gcf, 'Color', 'white')
plot([-t(2) 2*pi], t.*[1 1], 'Color', [1 1 1]*.6, 'LineWidth', 0.25); hold on
plot(t.*[1 1], [-t(2) 2*pi], 'Color', [1 1 1]*.6, 'LineWidth', 0.25)
mesh(taus{1}, taus{2}, taus{1}*0, 'EdgeColor', 'k', 'LineWidth', 1); hold on
plot(taus{1}, taus{2}, 'ko', 'MarkerFaceColor', 'w'); hold on
axis equal

Np = 8;
am1 = 20;
am2 = 3.5;
aa = gobjects(2);
cols = repmat([1 0 0], 2*Np+1, 1);
kk = 1;

xsh = pi+t(2)-(2*pi*3/Np+(ws(1)/ws(2))*(pi+t(2)));
for i=-Np:Np
    quiver(2*pi*i/Np+xsh, -t(2), ws(1)/vecnorm(ws)*am1, ws(2)/vecnorm(ws)*am1, '-', 'LineWidth', 1, 'Color', cols(kk,:));
    quiver(2*pi*i/Np+xsh, -t(2), ws(1)/vecnorm(ws)*am2, ws(2)/vecnorm(ws)*am2, '-', 'LineWidth', 1.5, 'Color', cols(kk,:));
    kk = kk+1;
end

tt = text(taus{1}(:)-0.4, taus{2}(:)-0.2, num2str((1:Nt^Nc)'), 'FontWeight', 'bold');

tics = (0:.25:1);
ticls = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
set(gca, 'XTick', tics*2*pi)
set(gca, 'XTickLabel', ticls)
set(gca, 'YTick', tics*2*pi)
set(gca, 'YTickLabel', ticls)

xlim([-t(2) 2*pi])
ylim([-t(2) 2*pi])
xlabel('Scaled Time $\tau_1$')
ylabel('Scaled Time $\tau_2$')

pti = 37;
mesh(taus{1}(pti)+[-1.5 -1.5;.5 .5]*t(2), taus{2}(pti)+[-1.5 .5;-1.5 .5]*t(2), zeros(2), 'FaceColor', 'none', 'LineStyle', '-.', 'EdgeColor', [0 1 0]*0.6, 'LineWidth', 2)

plot(taus{1}(pti)-t(2)*[0 1], taus{2}(pti)-t(2)*[1 0], '-', 'LineWidth', 1.5, 'Color', [0 1 0]*0.6)
yy = plot(taus{1}(pti)-ws(1)/sum(ws)*t(2), taus{2}(pti)-ws(2)/sum(ws)*t(2), 'kh', 'MarkerFaceColor', 'g');
mm = plot(taus{1}(pti)-t(2), taus{2}(pti), 'kh', 'MarkerFacecolor', 'b');
gg = plot(taus{1}(pti), taus{2}(pti)-t(2), 'kh', 'MarkerFacecolor', 'b');
rr = plot(taus{1}(pti), taus{2}(pti), 'kh', 'MarkerFaceColor', 'r') ;

if figsaves
    export_fig('./FIGS/E_qpgrid2.png', '-dpng')
end

xl = taus{1}(pti)+[-1.5 .5]*t(2);
yl = taus{2}(pti)+[-1.5 .5]*t(2);
xlim(xl); ylim(yl);
for i=1:length(tt)
    if ~any(i==[pti pti-1 pti-Nt pti-Nt-1])
        tt(i).Visible = 'off';
    else
        tt(i).Position = tt(i).Position+[0.3 0.25 0];
        tt(i).FontSize = 15;
    end
end
quiver(taus{1}(pti)-t(2), taus{2}(pti), t(2), 0, '-', 'LineWidth', 2, 'Color', [1 1 0]*0.5)
quiver(taus{1}(pti), taus{2}(pti)-t(2), 0, t(2), '-', 'LineWidth', 2, 'Color', [1 1 0]*0.5)
rr.MarkerSize = 15;
yy.MarkerSize = 15;
mm.MarkerSize = 15;
gg.MarkerSize = 15;

if figsaves
    export_fig('./FIGS/E_qpgrid2_zoom.png', '-dpng')
end

%%
[Nmat3, bpis, bpjs] = CONSTRUCTNMATHERE(ws, Nc, Nt, 3);
bpjs = cellfun(@(c) unique(c(:)), bpjs, 'UniformOutput', false);

figure(100)
clf()
set(gcf, 'Color', 'white')
plot([-t(2) 2*pi], t.*[1 1], 'Color', [1 1 1]*.6, 'LineWidth', 0.25); hold on
plot(t.*[1 1], [-t(2) 2*pi], 'Color', [1 1 1]*.6, 'LineWidth', 0.25)
mesh(taus{1}, taus{2}, taus{1}*0, 'EdgeColor', 'k', 'LineWidth', 1); hold on
plot(taus{1}, taus{2}, 'ko', 'MarkerFaceColor', 'w'); hold on
axis equal

Np = 8;
am1 = 20;
am2 = 3.5;
aa = gobjects(2);
cols = repmat([1 0 0], 2*Np+1, 1);
kk = 1;

xsh = pi+t(2)-(2*pi*3/Np+(ws(1)/ws(2))*(pi+t(2)));
for i=-Np:Np
    quiver(2*pi*i/Np+xsh, -t(2), ws(1)/vecnorm(ws)*am1, ws(2)/vecnorm(ws)*am1, '-', 'LineWidth', 1, 'Color', cols(kk,:));
    quiver(2*pi*i/Np+xsh, -t(2), ws(1)/vecnorm(ws)*am2, ws(2)/vecnorm(ws)*am2, '-', 'LineWidth', 1.5, 'Color', cols(kk,:));
    kk = kk+1;
end

tt = text(taus{1}(:)-0.4, taus{2}(:)-0.2, num2str((1:Nt^Nc)'), 'FontWeight', 'bold');

tics = (0:.25:1);
ticls = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
set(gca, 'XTick', tics*2*pi)
set(gca, 'XTickLabel', ticls)
set(gca, 'YTick', tics*2*pi)
set(gca, 'YTickLabel', ticls)

xlim([-t(2) 2*pi])
ylim([-t(2) 2*pi])
xlabel('Scaled Time $\tau_1$')
ylabel('Scaled Time $\tau_2$')

bi = 3;
ptis = [1:2 7];
plot([Nmat3(bpis{bi}(ptis), :)*taus{1}(:) taus{1}(bpis{bi}(ptis))]', [Nmat3(bpis{bi}(ptis), :)*taus{2}(:) taus{2}(bpis{bi}(ptis))]', 'b-.', 'LineWidth', 1.5);
gg = plot([Nmat3(bpis{bi}(ptis), :)*taus{1}(:) taus{1}(bpis{bi}(ptis))]', [Nmat3(bpis{bi}(ptis), :)*taus{2}(:) taus{2}(bpis{bi}(ptis))]', 'gh', 'MarkerFaceColor', 'g');
rr = plot(taus{1}(bpis{bi}(ptis)), taus{2}(bpis{bi}(ptis)), 'rh', 'MarkerFaceColor', 'r');

pti = bpis{bi}(ptis(1));
mesh(taus{1}(pti)+[-1.5 -1.5;1.5 1.5]*t(2), taus{2}(pti)+[-1.75 1.5;-1.75 1.5]*t(2), zeros(2), 'FaceColor', 'none', 'LineStyle', '-.', 'EdgeColor', [0 1 0]*0.6, 'LineWidth', 2)

if figsaves
    export_fig('./FIGS/E_qpgrid3.png', '-dpng')
end

xl = taus{1}(pti)+[-1.5 1.5]*t(2);
yl = taus{2}(pti)+[-1.75 1.5]*t(2);
xlim(xl); ylim(yl);
for i=1:length(tt)
    if ~any(i==[pti pti-1 pti-Nt pti-Nt-1])
        tt(i).Visible = 'off';
    else
        tt(i).Position = tt(i).Position+[0.3 0.25 0];
        tt(i).FontSize = 15;
    end
end
% quiver(taus{1}(pti)-ws(1)*deltau/vecnorm(ws), taus{2}(pti)-ws(2)*t(2)/vecnorm(ws), ws(1)*t(2)/vecnorm(ws), ws(2)*t(2)/vecnorm(ws), '-', 'LineWidth', 2, 'Color', [1 1 0]*0.5)
set(gg, 'MarkerSize', 15)
set(rr, 'MarkerSize', 15)

if figsaves
    export_fig('./FIGS/E_qpgrid3_zoom.png', '-dpng')
end

%%
figure(110)
clf()
set(gcf, 'Color', 'white')
plot([-t(2) 2*pi], t.*[1 1], 'Color', [1 1 1]*.6, 'LineWidth', 0.25); hold on
plot(t.*[1 1], [-t(2) 2*pi], 'Color', [1 1 1]*.6, 'LineWidth', 0.25)
mesh(taus{1}, taus{2}, taus{1}*0, 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', 'none'); hold on
plot(taus{1}, taus{2}, 'ko', 'MarkerFaceColor', 'w'); hold on
axis equal

Np = 8;
am1 = 20;
am2 = 3.5;
aa = gobjects(2);
cols = repmat([1 1 1]*0.6, 2*Np+1, 1);
kk = 1;

xsh = pi+t(2)-(2*pi*3/Np+(ws(1)/ws(2))*(pi+t(2)));
for i=-Np:Np
    quiver(2*pi*i/Np+xsh, -t(2), ws(1)/vecnorm(ws)*am1, ws(2)/vecnorm(ws)*am1, ...
        '-', 'LineWidth', 1, 'Color', cols(kk,:));
    quiver(2*pi*i/Np+xsh, -t(2), ws(1)/vecnorm(ws)*am2, ws(2)/vecnorm(ws)*am2, ...
        '-', 'LineWidth', 1.5, 'Color', cols(kk,:));
    kk = kk+1;
end

% tt = text(taus{1}(:)-0.4, taus{2}(:)-0.2, num2str((1:Nt^Nc)'), 'FontWeight', 'bold');

tics = (0:.25:1);
ticls = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
set(gca, 'XTick', tics*2*pi)
set(gca, 'XTickLabel', ticls)
set(gca, 'YTick', tics*2*pi)
set(gca, 'YTickLabel', ticls)

xlim([-t(2) 2*pi])
ylim([-t(2) 2*pi])
xlabel('Scaled Time $\tau_1$')
ylabel('Scaled Time $\tau_2$')

colos = DISTINGUISHABLE_COLORS(Nt);
for ti=1:Nt
    inds = bpis{ti};
    inds(1:(Nt-ti+1)) = inds((Nt-ti+1):-1:1);
    plot(taus{1}(inds), taus{2}(inds), '-', ...
        'MarkerFaceColor', colos(ti,:), 'Color', colos(ti,:), ...
        'LineWidth', 8)
%     if ti<Nt
%         mesh(taus{1}(ti:end,ti:end), taus{2}(ti:end,ti:end), taus{1}(ti:end,ti:end)*0+ti, 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', colos(ti,:)); hold on
%     end
end
% axis off
xlim(t([1 end])); ylim(t([1 end]))

if figsaves
    export_fig('./FIGS/E_qpgrid3_fronts.png', '-dpng')
end

%%
% Ncs = unique(fix(logspace(0, 5, 10)));  Ncs(1) = [];
% Nt = 4;
% 
% Ncss = [1 3 6 8];
% nnz3s = zeros(length(Ncss), 3);
% for ci=1:length(Ncss)
%     for ni=1:3
%         nnz3s(ci,ni) = nnz(CONSTRUCTNMATHERE(rand(1,Ncss(ci)), Ncss(ci), Nt, ni));
%     end
% end

Ncs = unique(fix(logspace(0, 5, 100)));  Ncs(1) = [];
Nt = 16;

figure(200)
clf()
loglog(Nt.^Ncs, (2./Nt).^Ncs, 'r-', 'LineWidth', 2); hold on
loglog(Nt.^Ncs, Ncs./(Nt.^Ncs), 'b-', 'LineWidth', 2); hold on
loglog(Nt.^Ncs, (2./Nt).^Ncs/2, 'g-.', 'LineWidth', 2); hold on
% loglog(Nt.^Ncss, nnz3s./(Nt.^(2*Ncss')), 'o-', 'LineWidth', 2)
xlabel('$N_t^{N_c}$')
ylabel('Sparsity ($\%$)')

legend('Approach 1', 'Approach 2', 'Approach 3')
grid on
% set(gca, 'XTick', Nt.^Ncss)
set(gcf, 'Color', 'white')
if figsaves
    export_fig('./FIGS/E_nmapsps.png', '-dpng')
end