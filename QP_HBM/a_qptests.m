clc
clear all
addpath('../ROUTINES/QUASIPERIODIC')

%%
Nc = 3;
Nt = 8;
Nhmax = fix(Nt/2)-1;

a0 = 10;
a1 = 1;
b1 = 3;
a2 = 2;
b2 = 4;

rng(1);
h1 = randi(Nhmax, 1, Nc);  % random h index
h2 = randi(Nhmax, 1, Nc);  % random h index

t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
taus = cell(Nc, 1);
[taus{:}] = ndgrid(t);

wtau = h1(1)*taus{1};
for ic=2:Nc
    wtau = wtau+h1(ic)*taus{ic};
end
y = a0 + a1*cos(wtau) + b1*sin(wtau);

wtau = h2(1)*taus{1};
for ic=2:Nc
    wtau = wtau+h2(ic)*taus{ic};
end
y = y + a2*cos(wtau) + b2*sin(wtau);

%% FFTn
iof = fix(Nt/2)+1;
i0 = num2cell(repmat(iof, 1, Nc));
i1 = num2cell(iof+h1);
j1 = num2cell(iof-h1);

yf = fftshift(fftn(y))*2/(Nt^Nc);
yf(i0{:}) = yf(i0{:})/2;  % zeroth harmonic

yf(i1{:})
yf(j1{:})

%%
Y = QPFOURIERCOEFF(y(:), [zeros(1, Nc); h1; h2]);
yt = QPTIMETRANS(Y, [zeros(1, Nc); h1; h2], Nt);
max(abs(yt(:)-y(:)))

%%
ws = rand(1,Nc);
% [E, dEdw] = QPHARMONICSTIFFNESS(1,1,1,ws, [zeros(1,Nc);h1;h2])
QPHARMONICSTIFFNESS(1,1,1,ws, [zeros(1,Nc);h1;h2])

%%
taue = rand(4,Nc);
Ns = QPTIMEINTERP(taue, [zeros(1,Nc);h1; h2])

%%
h = [zeros(1, Nc); h1; h2];
Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
DFTm = QPTIMETRANS(eye(Nhc), h, Nt);