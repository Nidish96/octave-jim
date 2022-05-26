clc
clear all
addpath('../../ROUTINES/')
addpath('../../ROUTINES/FEM/BEAMS')

%%
geom.E = 2e11;
geom.a = 0.013;
geom.b = 0.001;
geom.A = geom.a*geom.b;
geom.I = geom.a*geom.b^3/12;
geom.rho = 7800;
geom.L = 0.3;

Ne = 100;
[Me, Ke] = EBBEAM_MATS(geom.rho, geom.E, geom.A, geom.I, geom.L/Ne);

M = zeros((Ne+1)*3);
K = zeros((Ne+1)*3);
for e=1:Ne
    M((e-1)*3+1:(e+1)*3, (e-1)*3+1:(e+1)*3) = M((e-1)*3+1:(e+1)*3, (e-1)*3+1:(e+1)*3) + Me;
    K((e-1)*3+1:(e+1)*3, (e-1)*3+1:(e+1)*3) = K((e-1)*3+1:(e+1)*3, (e-1)*3+1:(e+1)*3) + Ke;
end
Lb = eye((Ne+1)*3);
Lb(:, [1:3 end-2:end]) = [];

Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;

%%
[Vb, Db] = eig(Kb, Mb);
Wb = sqrt(diag(Db))/2/pi;

mi = 1;
di = 2;

fint = K*(Lb*Vb);

figure(1)
clf()
plot(Lb(di:3:end,:)*Vb(:, mi), '.-')
yyaxis right
plot(K(di:3:end, :)*(Lb*Vb(:,mi)), '.-')
% plot(Lb(di:3:end,:)*(Kb*Vb(:,mi)), 'o-')