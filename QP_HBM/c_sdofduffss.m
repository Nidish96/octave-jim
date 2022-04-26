clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/QUASIPERIODIC')

anim = false; 

%%
m = 1;
c = 0.5;
k = 4;
bt = 0.15;   % 0.5
fnl = @(t,u,ud) deal(bt*u.^3, 3*bt*u.^2, zeros(size(u)));
% fnl = @(t,u,ud) deal(bt*u.^3 + c*ud, 3*bt*u.^2, c*ones(size(u)));

Nc = 2;  % Number of components
Nhmax = 10;  % Number of harmonics
%% Harmonic Selection
hall = cell(1, Nc);
[hall{:}] = ndgrid(-Nhmax:Nhmax);
hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
h = hall(sum(abs(hall),2)<=Nhmax & sum(hall,2)>=0,:);
% h = hall(abs(sum(hall,2))<=Nhmax & sum(hall,2)>=0,:);

h(sum(h,2)==0 & h(:,1)<=0, :) = [];
h = [zeros(1,Nc); h];

figure(1)
clf()
plot(hall(:,1), hall(:,2), 'ko', 'MarkerFaceColor', 'w'); hold on
plot(h(:,1), h(:,2), '*'); hold on
grid on
axis equal
legend('All Harmonics', 'Selected Harmonic', 'Location', 'northoutside')

%% Setup Model
GM = MDOFGEN(m, k, c, 1.0);
GM = GM.SETNLFUN(1+3, 1.0, fnl);

%% Forcing
% ws = [pi sqrt(2)];
ws = [sqrt(k/m) pi.^(1:Nc-1)];
% rng(1);
% hid = randi(length(h)-1,2,1)
% hid = [1; 2];

hid = [find(h(:,1)==0 & h(:,2)==1); find(h(:,1)==1 & h(:,2)==0)]-1;
hfrc = h(1+hid, :);

amps = 20*ones(size(hid));  % 20
fext = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

%% Transient simulation
T0 = 0;
T1 = 400;

fsamp = 60;

opts = struct('Display', 'waitbar');
[T, U, Ud, Udd] = GM.HHTAMARCH(T0, T1, 1/fsamp, 0, 0, fext, opts);

%%
Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));

E = QPHARMONICSTIFFNESS(GM.M, GM.C, GM.K, ws, h);
Fl = zeros(Nhc, 1);
Fl(1+(hid-1)*2+1) = amps;

D1 = QPHARMONICSTIFFNESS(0, 1, 0, ws, h);  % Time derivative matrix

X0 = E\Fl;
Nt = 64;

% [R0, dR0] = GM.QPHBRESFUN([X0; 1], ws, Fl, h, Nt, eps);
% %%

% opt = struct('Display', true, 'ITMAX', 200);
% X = NSOLVE(@(U) GM.QPHBRESFUN([U; 1], ws, Fl, h, Nt, eps), X0, opt);

fopt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
X = fsolve(@(U) GM.QPHBRESFUN([U; 1], ws, Fl, h, Nt, eps), X0, fopt);

Xd = D1*X;
Xdd = D1^2*X;
Ns = QPTIMEINTERP(T(:).*ws, h);
%%
figure(2)
clf()
plot(T, U, 'b-', 'LineWidth', 2)
grid on; hold on
% plot(T, Y, '.-');
plot(T, Ns*X, 'r-', 'LineWidth', 2)

figure(3)
clf()
plot(U, Ud, 'b-', 'LineWidth', 2); grid on; hold on
plot(Ns*X, Ns*Xd, 'r-', 'LineWidth', 1)

%% Torus
if Nc==2
    Nt = 64;
    x   = reshape(QPTIMETRANS(X, h, Nt), repmat(Nt, 1, Nc));
    xd  = reshape(QPTIMETRANS(D1*X, h, Nt), repmat(Nt, 1, Nc));
    xdd = reshape(QPTIMETRANS(D1^2*X, h, Nt), repmat(Nt, 1, Nc));
    
    x   = [x   x(:,1)  ; x(1,:)   x(1,1)];
    xd  = [xd  xd(:,1) ; xd(1,:)  xd(1,1)];
    xdd = [xdd xdd(:,1); xdd(1,:) xdd(1,1)];
    figure(4)
    clf()
    % plot3(x, xd, xdd, '-', 'LineWidth', 1)
    surf(x, xd, xdd, 'EdgeColor', 'r'); hold on
    plot3(U(1000:end), Ud(1000:end), Udd(1000:end), 'b-')
    grid on
end

%% Frequency Content
wndw = blackmanharris(length(T));

[frqs, uf] = FFTFUN(T', U'.*wndw);
[~, xf] = FFTFUN(T', (Ns*X).*wndw);

figure(5)
clf()
semilogy(frqs, abs(uf), 'b-'); hold on
% plot(frqs, abs(xf), 'r-')
stem(abs(h*ws(:)/2/pi), abs([X(1); X(2:2:end)+1j*X(3:2:end)]), 'filled', 'k')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 5])

%% Animate
if Nc==2 && anim
    Nt = 64;
    x   = reshape(QPTIMETRANS(X, h, Nt), repmat(Nt, 1, Nc));
    xd  = reshape(QPTIMETRANS(D1*X, h, Nt), repmat(Nt, 1, Nc));
    xdd = reshape(QPTIMETRANS(D1^2*X, h, Nt), repmat(Nt, 1, Nc));
    
    x   = [x   x(:,1)  ; x(1,:)   x(1,1)];
    xd  = [xd  xd(:,1) ; xd(1,:)  xd(1,1)];
    xdd = [xdd xdd(:,1); xdd(1,:) xdd(1,1)];
    
    figure(4)
    clf()
    surf(x, xd, xdd, 'EdgeColor', 'r'); hold on
	plot3(U(1000:end), Ud(1000:end), Udd(1000:end), 'b-')
	grid on    
    xl = xlim; yl = ylim; 
    
    for ti=1:100:length(T)
        clf()
        % plot3(x, xd, xdd, '-', 'LineWidth', 1)
        surf(x, xd, xdd, 'EdgeColor', 'none', 'FaceAlpha', 1); hold on
        plot3(U(1:ti), Ud(1:ti), Udd(1:ti), '-', 'Color', [1 1 1]*0.6)
        plot3(U(ti), Ud(ti), Udd(ti), 'bo', 'MarkerFaceColor', 'b')
        grid on
        xlim(xl); ylim(yl);
        title(sprintf('Frame %d/%d', ti, length(T)))
        pause(0.01)
    end
end