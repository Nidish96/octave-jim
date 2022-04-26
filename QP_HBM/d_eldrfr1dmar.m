clc
clear all
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

%%
Nt = 128;

t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
u = sin(t);

kt = 1;
muN = 0.5;

%% Scheme 1
f = kt*u;
fprev = f(end);

figure(1)
clf()
aa = gobjects(11,1);
aa(1) = plot(t-2*pi, u, '-', 'Color', [1 1 1]*0.6, 'LineWidth', 1.5); hold on
aa(2) = plot(t-2*pi, f, '-', 'LineWidth', 2);
legend(aa(1), 'Lin. Ref.')
legend(aa(2), 'Iter. 0')
grid on
tic
for it=1:10
    plot(t + (it-1)*2*pi, u, '-', 'Color', [1 1 1]*0.6, 'LineWidth', 1.5); hold on
    
    for ti=1:Nt
        tmi = mod(ti-1-1, Nt)+1;
        
        fsp = kt*(u(ti)-u(tmi)) + f(tmi);
        if abs(fsp)<muN
            f(ti) = fsp;
        else
            f(ti) = muN*sign(fsp);
        end
    end
    aa(it+2) = plot(t + (it-1)*2*pi, f, '-', 'LineWidth', 2); hold on
    legend(aa(it+2), sprintf('Iter. %d', it))
    title(sprintf('Iter %d: %e', it, abs(fprev-f(end))))
    if abs(fprev-f(end))<eps
        break
    end
    fprev = f(end);
end
toc
aa(it+3) = plot([xlim nan xlim], [1 1 nan -1 -1]*muN, '-', 'Color', [1 1 1]*0.6);
legend(aa(it+3), 'Slip limit')
legend(aa(1:it+3), 'Location', 'southwest')
xlabel('Scaled Time $(\tau=t/\omega)$')
ylabel('Force (N)')
set(gcf, 'Color', 'white')
title('')
tics = (-2:.5:it)*2;
set(gca, 'XTick', tics*pi)
set(gca, 'XTickLabel', [num2str(tics(:)) repmat('\pi', length(tics), 1)])
export_fig('./FIGS/1DMARCH.png', '-dpng')
%% Scheme 2 (half step interpolation): faster convergence when closer to tmi
wgt = 0.05;

figure(1)
clf()

f = kt*u;
fprev = f(end);

plot(t-2*pi, u, 'Color', [1 1 1]*0.6); hold on
plot(t-2*pi, f, 'b-')
grid on
errs = [];
for it=1:1000
    plot(t + (it-1)*2*pi, u, 'Color', [1 1 1]*0.6); hold on
    
    for ti=1:Nt
        tmi = mod(ti-1-1, Nt)+1;
        Nsf = [1-wgt wgt];
        
        fsp = kt*(u(ti) - Nsf*u([tmi ti])) + Nsf*f([tmi ti]);
        if abs(fsp)<muN
            f(ti) = fsp;
        else
            f(ti) = muN*sign(fsp);
        end
    end
    
    plot(t + (it-1)*2*pi, f, 'b-'); hold on
    errs = [errs; abs(fprev-f(end))];
    title(sprintf('Iter %d: %e', it, errs(end)))
    xlim([it-2 it]*2*pi);
    if errs(end)<eps
        break
    end
    fprev = f(end);
    
end

figure(3)
semilogy(errs+eps/2, 'o-')

%% Scheme 3 (half step extrapolation not involving present point): very good 
wgt = 1.9;

figure(1)
clf()

f = kt*u;
fprev = f(end);

plot(t-2*pi, u, 'Color', [1 1 1]*0.6); hold on
plot(t-2*pi, f, 'b-')
grid on
errs = [];
for it=1:1000
    plot(t + (it-1)*2*pi, u, 'Color', [1 1 1]*0.6); hold on
    
    for ti=1:Nt
        tmi1 = mod(ti-1-1, Nt)+1;
        tmi2 = mod(ti-2-1, Nt)+1;
        
        Nsf = [1-wgt wgt];
        
        fsp = kt*(u(ti) - Nsf*u([tmi2 tmi1])) + Nsf*f([tmi2 tmi1]);
        if abs(fsp)<muN
            f(ti) = fsp;
        else
            f(ti) = muN*sign(fsp);
        end
    end
    
    plot(t + (it-1)*2*pi, f, 'b-'); hold on
    errs = [errs; abs(fprev-f(end))];
    title(sprintf('Iter %d: %e', it, errs(end)))
    xlim([it-2 it]*2*pi);
    if errs(end)<eps
        break
    end
    fprev = f(end);
    
end

figure(3)
semilogy(errs+eps/2, 'o-')

%% Scheme 4 (second order formula)
wgt = 1.5;

figure(1)
clf()

f = kt*u;
fprev = f(end);

plot(t-2*pi, u, 'Color', [1 1 1]*0.6); hold on
plot(t-2*pi, f, 'b-')
grid on
errs = [];
for it=1:1000
    plot(t + (it-1)*2*pi, u, 'Color', [1 1 1]*0.6); hold on
    
    for ti=1:Nt
        tpi1 = mod(ti+1-1, Nt)+1;
        tmi1 = mod(ti-1-1, Nt)+1;
        tmi2 = mod(ti-2-1, Nt)+1;
        
        Nsf = [1-wgt wgt];
        
        fsp = kt*(u(ti) - u(tmi2)) + f(tmi2);
        if abs(fsp)<muN
            f(ti) = fsp;
        else
            f(ti) = muN*sign(fsp);
        end
    end
    
    plot(t + (it-1)*2*pi, f, 'b-'); hold on
    errs = [errs; abs(fprev-f(end))];
    title(sprintf('Iter %d: %e', it, errs(end)))
    xlim([it-2 it]*2*pi);
    if errs(end)<eps
        break
    end
    fprev = f(end);
    
end

figure(3)
semilogy(errs+eps/2, 'o-')

%% Scheme 5 (implicit march): works but has phase inaccuracies for large wgt
wgt = 0.5;

figure(1)
clf()

f = kt*u;
fprev = f(end);

plot(t-2*pi, u, 'Color', [1 1 1]*0.6); hold on
plot(t-2*pi, f, 'b-')
grid on
errs = [];

opts = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
tic
for it=1:1000
    plot(t + (it-1)*2*pi, u, 'Color', [1 1 1]*0.6); hold on
    
    for ti=1:Nt
        tmi = mod(ti-1-1, Nt)+1;
        Nsf = [1-wgt wgt];
        
        
        fun = @(fr) deal(fr-(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])*(abs(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])<muN) - muN*sign(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])*(abs(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])>=muN), ...
            1-Nsf(2).*(abs(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])<muN));
        
        f(ti) = fsolve(fun, f(ti), opts);
    end
    
    plot(t + (it-1)*2*pi, f, 'b-'); hold on
    errs = [errs; abs(fprev-f(end))];
    title(sprintf('Iter %d: %e', it, errs(end)))
    xlim([it-2 it]*2*pi);
    if errs(end)<eps
        break
    end
    fprev = f(end);
    
end
toc
figure(3)
semilogy(errs+eps/2, 'o-')
fsol = f;

%% Scheme 6 (one-shot implicit march): works
wgt = 0.5;

figure(1)
clf()

f = kt*u;
fprev = f(end);

plot(t-2*pi, u, 'Color', [1 1 1]*0.6); hold on
plot(t-2*pi, f, 'b-')
grid on
errs = [];

opts = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
optsr = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', false);
tic
% 1-shot march
Dm = diag(ones(Nt-1,1)*(1-wgt),-1) + eye(Nt)*wgt; Dm(1, end) = 1-wgt;

fspr = @(fr) fr - (kt*(eye(Nt)-Dm)*u+Dm*fr).*(abs(kt*(eye(Nt)-Dm)*u+Dm*fr)<muN) - muN*sign(kt*(eye(Nt)-Dm)*u+Dm*fr).*(abs(kt*(eye(Nt)-Dm)*u+Dm*fr)>=muN);
fsp = @(fr) deal(fr - (kt*(eye(Nt)-Dm)*u+Dm*fr).*(abs(kt*(eye(Nt)-Dm)*u+Dm*fr)<muN) - muN*sign(kt*(eye(Nt)-Dm)*u+Dm*fr).*(abs(kt*(eye(Nt)-Dm)*u+Dm*fr)>=muN), ...
    eye(Nt) - Dm.*(abs(kt*(eye(Nt)-Dm)*u+Dm*fr)<muN));

fs = fsolve(fsp, f, opts);
% fsr = fsolve(fspr, f, optsr);
toc

tic
for it=1:1000
%     plot(t + (it-1)*2*pi, u, 'Color', [1 1 1]*0.6); hold on
    
    for ti=1:Nt
        tmi = mod(ti-1-1, Nt)+1;
        Nsf = [1-wgt wgt];
        
        
        fun = @(fr) deal(fr-(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])*(abs(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])<muN) - muN*sign(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])*(abs(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])>=muN), ...
            1-Nsf(2).*(abs(kt*(u(ti)-Nsf*u([tmi ti]))+Nsf*[f(tmi);fr])<muN));
        
        f(ti) = fsolve(fun, f(ti), opts);
    end
    
%     plot(t + (it-1)*2*pi, f, 'b-'); hold on
    errs = [errs; abs(fprev-f(end))];
%     title(sprintf('Iter %d: %e', it, errs(end)))
%     xlim([it-2 it]*2*pi);
    if errs(end)<eps
        break
    end
    fprev = f(end);
    
end
toc
% figure(3)
% semilogy(errs+eps/2, 'o-')
fsol = f;

figure(1)
clf()
plot(t, fs, '.-'); hold on; grid on; 
plot(t, fsol, 'o-'); 
plot(t, u, '.-'); 
plot(t(Nt/4+[1 1]), ylim, '-'); 
plot(t(3*Nt/4+[1 1]), ylim, '-'); hold off
xlabel('Time')

%% Scheme 7 - derivative solve - not working
Dm = real(ifft(1j*[0:(Nt/2-1) (-Nt/2):-1]'.*fft(eye(Nt))));  % Fourier Differentiation Matrix
% Dm = real(ifft([0; 1./(1j*[1:(Nt/2-1) (-Nt/2):-1]')].*fft(eye(Nt))));  % Fourier Integration Matrix

Dm = spdiags(ones(Nt,1).*[-1 1], [-1 1], Nt, Nt); Dm(1, end) = -1; Dm = Dm*Nt/4/pi;

Dmi = pinv(full(Dm));

u = sin(t);
ud = cos(t);

resfun = @(f) deal(f-Dmi*((kt*Dm*u).*(abs(f)<muN)), eye(Nt));
% resfun = @(f) f-Dm*((kt*ud).*(abs(f)<muN));

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
fsol = fsolve(resfun, kt*u, opt);

opt = struct('Display', true);
fsol = NSOLVE(resfun, kt*u, opt);