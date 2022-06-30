clc
clear all

%% Checking out Lagrangian Interpolation
ptsb = (cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^Nc-1))))*2-1);  % using binary for the construction of points on a unit square
xis = ptsb*0.1;

fun = @(x,y) exp(-(x.^2+y.^2)).*x;

u = fun(ptsb(:, 1), ptsb(:, 2));

Lm = prod(diff(xis([1 end],:)));  % Area of unit cube

xyq = [0.01 0.01];
Nsf = prod(abs(xis(end:-1:1,:)-xyq), 2)'/Lm
uhatq =  Nsf*fun(ptsb(:,1), ptsb(:,2))
fq = fun(xyq(1), xyq(2))

%%
N = 10;
[x, y] = meshgrid(linspace(-2, 2, N), linspace(-2, 2, N));

[xx, yy] = meshgrid(linspace(min(xis(:,1)), max(xis(:,1)), N)*5, 5*linspace(min(xis(:,2)), max(xis(:,2)), N));
Nsx = squeeze(prod(xis(end:-1:1,:)-permute([xx(:), yy(:)], [3 2 1]),2))'/Lm;

figure(1)
clf()
surf(xx, yy, fun(xx, yy)); hold on
plot3(xx(:), yy(:), Nsx*fun(xis(:,2), xis(:,1)), 'ko', 'MarkerFaceColor', 'k')