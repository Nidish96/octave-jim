% clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/QUASIPERIODIC')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

anim = false; 
savdat = true;
%% Parameter Sets for Examples
parmslist = [struct('m', 1, 'k', 1, 'c', 1e-3, 'kt', 3, 'muN', 1.5, ...
    'ws', @(Nc) sqrt(2:Nc+1), 'amps', 0.375, 'epN', 1e-1, 'ss', 1.0);  % Example 1
    struct('m', 1, 'k', 1, 'c', 1e-3, 'kt', 3, 'muN', 1.5, ...
    'ws', @(Nc) sqrt(2:Nc+1), 'amps', 5, 'epN', 1e-1, 'ss', 1.0);  % Example 1
    struct('m', 1, 'k', 80, 'c', 2*0.1e-2*sqrt(80), 'kt', 100, 'muN', 10, ...
    'ws', @(Nc) sqrt(((1:Nc)-1)*100+2*10.^(1:Nc)), 'amps', 10, 'epN', 1e-1, 'ss', 2.0);  % Example 2
    struct('m', 1, 'c', 0.01, 'k', 1, 'kt', 3, 'muN', 0.5, ...
    'ws', @(Nc) sqrt(2:Nc+1), 'amps', 5, 'epN', 1e-6, 'ss', 1.0);    
    ];

%%
ex = 1;
Nc = 2;
ws = parmslist(ex).ws(Nc);

MDL = MDOFGEN(parmslist(ex).m, parmslist(ex).k, parmslist(ex).c, 1.0);
MDL = MDL.SETNLFUN(2+3, 1, @(t, u, varargin) JENKNL(t, u, parmslist(ex).kt, parmslist(ex).muN, varargin{:}), [], 4);

%% Harmonic Selection
Nhmax = 1;
Nt = 128;
h = HSEL(Nhmax, ws);
Nhc = sum(all(h==0,2) + 2*any(h~=0,2));

[rinds0,zinds,hinds,rinds,iinds] = HINDS(1, h);

%% Excitation
Fl = zeros(Nhc,1);
Fl(rinds(1:2)) = parmslist(ex).amps;

U0 = QPHARMONICSTIFFNESS(MDL.M, MDL.C, MDL.K, ws, h)\Fl;

JNL = zeros(Nhc, Nhc+Nc);
[FNL, JNL(:,1:Nhc), JNL(:,Nhc+1:end)] = MDL.QPNLEVAL([U0;ws(:)], h, Nt, 1e-10);

JNUM = zeros(Nhc, Nhc+Nc);
parfor hi = 1:Nhc+Nc
    hv = zeros(Nhc+Nc,1);
    hm = 1e-6;
    hv(hi) = 1;
    FNLp = MDL.QPNLEVAL([U0;ws(:)]+hv*hm, h, Nt, 1e-10);
    FNLm = MDL.QPNLEVAL([U0;ws(:)]-hv*hm, h, Nt, 1e-10);
    hv(hi) = 0;

    JNUM(:, hi) = (FNLp-FNLm)/(2*hm);
    fprintf('%d/%d\n', hi, Nhc+Nc);
end