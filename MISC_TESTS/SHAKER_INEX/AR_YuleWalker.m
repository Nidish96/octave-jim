clc
clear all

%% Understanding autocorr
fsamp = 2^7;
T0 = 0;
T1 = 1;

t = (T0:1/fsamp:T1)';
Nt = length(t);
x = cos(2*pi*5*t);
x = rand(Nt,1);

Rxx = zeros(Nt-1, 1);
for i=1:Nt-1
    Rxx(i) = mean(x(1:end-(i-1)).*x(i:end));
end

plot([Rxx(end:-1:2); Rxx], '.-')