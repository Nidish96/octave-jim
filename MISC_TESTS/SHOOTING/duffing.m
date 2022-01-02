clc
clear all

%% 
m = 1.0;
k = 4.0;
bt = 2;
c = 2*0.01*sqrt(k/m)*m;

%%
MDL = MDOFGEN(m, k, c, 1);

fnl = @(t, u, ud) deal(bt*u.^3, 3*bt*u.^2, zeros(size(u)));
MDL = MDL.SETNLFUN(1+3, 1.0, fnl);

%%
Wfrc = 2.0;
f0 = 1;
Fex = @(t) f0*cos(Wfrc*t);

Ncyc = 1;
Ntcyc = 128;
Tmax = 2*pi/Wfrc*Ncyc;

opts = struct('Display', 'waitbar');
[T, U, Ud, ~, ~, PHI] = MDL.HHTAMARCH(0, Tmax, Tmax/(Ntcyc*Ncyc), ...
    0, 0, Fex, opts);
%% Single Frequency Solve
opt = struct('Display', true);

famp = 10;
xs = NSOLVE(@(u) MDL.SHOOTINGRES([u; Wfrc], famp, 1, Ntcyc, Ncyc), [0; 0], opt);

% %%
[~, ~, T, U, Ud] = MDL.SHOOTINGRES([xs; Wfrc], famp, 1, Ntcyc, Ncyc);

% plot(T, U, '.-')

%% Continuation
Copt = struct('Nmax', 1000, 'DynDscale', 1, ...
    'arclengthparm', 'orthogonal', 'angopt', 5e-2, ...
    'parmjac', false);

famp = 0.02;
Wst = 1.5;
Wen = 2.5;
dw = 0.1;

resfun = @(uw) MDL.SHOOTINGRES(uw, famp, 1, Ntcyc, Ncyc);
sol = CONTINUE(resfun, [0; 0], Wst, Wen, dw, Copt);

%%
Amps = zeros(size(sol,2), 1);
stab = zeros(size(sol,2), 1);
eVs = zeros(size(sol,2), 2);

figure(3)
clf()

for iw=1:size(sol,2)
    [~, ~, T, U, Ud, PHI] = MDL.SHOOTINGRES(sol(:, iw), famp, 1, Ntcyc, Ncyc);
    Amps(iw) = max(U);
    
    eVs(iw, :) = eig(PHI);
    stab(iw) = any(abs(eig(PHI))>1);
    
    fprintf('%d/%d\n', iw, size(sol,2))
    
    figure(3)
    clf()
    plot(T*sol(end,iw), U, '.-'); hold on
    if stab(iw)==0
        title('Stable')
    else
        title('Unstable')
    end
    grid on
    ylim(0.5*[-1 1])
    pause(0.1)
    
end
%%

figure(1)
% hold on
clf()
plot(sol(end,:)./(1-stab'), Amps./(1-stab'), '-'); hold on
plot(sol(end, :)./stab', Amps./stab', '--')
xlabel('Frequency (rad/s)')
ylabel('Max Amplitude')


%%
figure(2)
clf()

plot(eVs(:, 1), '.-'); hold on
plot(eVs(:, 2), '.-'); 

t = linspace(0, 2*pi, 100);
plot(cos(t), sin(t), 'k-')

xlabel('Real')
ylabel('Imag')

grid on
axis equal