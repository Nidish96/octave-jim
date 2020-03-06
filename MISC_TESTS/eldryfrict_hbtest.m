clc
clear all
addpath('../ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/HARMONIC')

%% Parameters
Np = 1;  % Number of points
kxynmu = repmat([1 1 2 0.5], 1, Np)';

h = [0 1];  Nhc = sum(h==0)+2*sum(h~=0);
Nt = 128;

rand("seed", 1);
UXYN = rand(3*Np*Nhc,1);
N0 = zeros(Np,1)+1;

[FXYN, JXYN] = ELDRYFRICT_HB(UXYN, h, Nt, kxynmu, N0);

hm = 1e-4;
hv = zeros(size(UXYN));

JXYN_num = zeros(size(JXYN));
for i=1:length(UXYN)
  hv(i) = 1;
  FXYNp = ELDRYFRICT_HB(UXYN+hv*hm, h, Nt, kxynmu, N0);
  FXYNm = ELDRYFRICT_HB(UXYN-hv*hm, h, Nt, kxynmu, N0);

  JXYN_num(:, i) = (FXYNp-FXYNm)/(2*hm);
  hv(i) = 0;
  fprintf('Done %d/%d\n', i, length(UXYN));
end

disp(max(max(abs(JXYN-JXYN_num))))
