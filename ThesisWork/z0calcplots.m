clc
clear all
addpath('../ROUTINES/export_fig')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

%%
E = 2e11; nu = 0.30;
Estar = E/(1-nu^2);
Gstar = E/((1+nu)*(2-nu));

Rad = 0.1e-3;
lam = 1e6;
Nasps = 1000;

cn = Estar*sqrt(pi*Rad./lam.^3)*Nasps;

us = linspace(-1, 1, 1000)*(25/lam);

z0 = log(Nasps)/lam;
funcu = @(u) 4/3*Estar*sqrt(Rad/lam^3)*(3/4*sqrt(pi)*exp(lam*(z0+u)).*erf(abs(sqrt(lam*(z0+u))))-3/2*sqrt(lam*(z0+u))).*(u>=-z0);
funch = @(u) 4/3*Estar*sqrt(Rad)*(z0+u).^1.5;

figure(1)
clf()
set(gcf, 'Color', 'white')
a3 = semilogy(us*1e6, real(funch(us)), 'r-', 'LineWidth', 2); hold on;
a1 = semilogy(us*1e6, cn*exp(lam*us), 'b-', 'LineWidth', 2); 
a4 = semilogy(us*1e6, real(funcu(us)), 'g-.', 'LineWidth', 4); 
a2 = plot(-z0*[1 1]*1e6, ylim, 'k-', 'LineWidth', 2); 

legend([a1 a2 a3 a4], 'Analytical Model', 'Height $z_0$', ...
       'Single Asperity Model', 'Lu \& O''Shea (2012)', 'Location', 'southeast'); 

grid on;
xlabel('Normal Displacement ($\mu m$)')
ylabel('Normal Force (N)')

export_fig('./REPS/FIGS/z0modif.png', '-dpng')