clc
clear all
addpath('../ROUTINES/QUADRATURE/')

%% Gauss-Lagrange Quadrature (accuracy up to 2n-1)
% Represents:
% \int_{-1}^1 f(x) dx = \sum_i wi f(xi)
Nq = 5;
[xi, wi] = LGWT(Nq);

%% Gauss-Lobatto Quadrature (end points included; accuracy is up to 2n-3)
% Represents:
% \int_{-1}^1 f(x) dx = \sum_i wi f(xi)
Nq = 5;
[xi, wi] = LGLWT(Nq);

%% Probabilists' Gauss-Hermite Quadrature (accurate up to 2n-1)
% Represents:
% \int_{\inf}^{\inf} f(x) \frac{1}{\sqrt{2\pi}}exp(-x^2/2) dx = \sum_i wi f(xi)

Nq = 5;
[xi, wi] = GPHWT(Nq);

%% Gauss-Laguerre Quadrature (accurate up to 2n-1)
% Represents:
% \int_0^{\inf} f(x) exp(-x) dx = \sum_i wi f(xi)

Nq = 5;
[xi, wi] = LAGWT(Nq);

