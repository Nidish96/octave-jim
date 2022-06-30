% clc
clear all
addpath('../../ROUTINES/SOLVERS/')

%%
sc = 100;

Copt = struct('Nmax', 200, 'Display', 1, 'DynDscale', 1, 'solverchoice', 1, 'angopt', 1e-1);
Copt.Dscale = [1.0; sc];

xstart = -1*sc;
xend = 1*sc;
dx = 0.15;

Xs = CONTINUE(@(yx) resfun(yx, sc), -2, xstart, xend, dx, Copt);

figure(1)
clf()
% for i=1:size(Xs,2)
%     clf()
%     plot(Xs(2,1:i), Xs(1,1:i), '.-', Xs(2,1), Xs(1,1), 'ko'); 
%     xlim([xstart xend])
%     ylim([-1 1]*1.5)
%     pause(0.1)
% end
plot(Xs(2,:), Xs(1,:), '.-', Xs(2,1), Xs(1,1), 'ko'); 
xlim([xstart xend])
ylim([-1 1]*1.5)

% opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% fsolve(@(y) resfun([y; xstart]), -2, opts)

function [res, dresdy, dresdx] = resfun(yx, sc)
    x = yx(2);
    y = yx(1);

    res = x-sc*(y^3-y);
    dresdy = -sc*(3*y^2-1);
    dresdx = 1.0;
end