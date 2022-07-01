% QP marching approach for force estimation
clc
clear all
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

plotfigs = true;

Nmtype = 3;  % 1-interpolation ; 2-fdm; 3-exclusive interpolation
%% Frequency configuration
Nt = 16;
Nc = 2;
Nhmax = 5;

ws = sqrt(1:Nc);
% ws = [pi 1.0];

% ws(1) = ws(1)*0.01;

%% Harmonic Selection
hall = cell(1, Nc);
[hall{:}] = ndgrid(-Nhmax:Nhmax);
hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
h = hall(abs(sum(hall,2))<=Nhmax & sum(hall,2)>=0,:);

h(sum(h,2)==0 & h(:,1)<=0, :) = [];
h = [zeros(1,Nc); h];

%%
t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
deltau = t(2);
taus = cell(Nc, 1);
[taus{:}] = ndgrid(t);

tausn = cell(Nc, 1);
[tausn{:}] = ndgrid(1:Nt);

%% Scheme 1: Build Interpolation Shape Functions
dt_vec = ws*deltau/vecnorm(ws);  % vector corresponding to deltatau amplitude in real time dxn

ptsb = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^Nc-1))));  % using binary for the construction of points on a unit square
oppi = (2^Nc):-1:1;  % diagonally opposite points are retrieved using the binary inverses
xis = ptsb*(-deltau);  % coordinates in "tau" space relative to origin

Lm = deltau^Nc;  % Lebesgue Measure of each cell in tau space
Nsf = prod(abs(xis(oppi,:)-(-dt_vec)), 2)'/Lm;  % Lagrange Shape Functions to interpolate previous point in the cell diagonally behind (in tau space)
    
ijs = cell2mat(cellfun(@(c) c(:), tausn, 'UniformOutput', false)');  % indices of all points
evijs = mod(repmat(ijs, 1, 1, 2^Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the cell that is diagonally behind
evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

% Build sparse interpolation matrix
tic
Nmat1 = sparse(repmat((1:Nt^Nc)', 1, 2^Nc), evns, repmat(Nsf, Nt^Nc, 1));
toc
% note to self: Use Nmat as "Dm" is used in scheme 6 in "d_eldrfr1dmar.m"

%% Scheme 2: Build FD shape functions
ptsb = eye(Nc);

Nsf = ws/sum(ws);
    
ijs = cell2mat(cellfun(@(c) c(:), tausn, 'UniformOutput', false)');  % indices of all points
evijs = mod(repmat(ijs, 1, 1, Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the cell that is diagonally behind
evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

% Build sparse interpolation matrix
tic
Nmat2 = sparse(repmat((1:Nt^Nc)', 1, Nc), evns, repmat(Nsf, Nt^Nc, 1));
toc
% note to self: Use Nmat as "Dm" is used in scheme 6 in "d_eldrfr1dmar.m"

%% Scheme 3: Use the "previous plane" so as to make the interpolation explicit
ptsb = ones(1, Nc);  % diagonally opposite point

ijs = cell2mat(cellfun(@(c) c(:), tausn, 'UniformOutput', false)');  % indices of all points
evijs = mod(repmat(ijs, 1, 1, 1) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the cell that is diagonally behind

evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

% Build interpolant matrix
dt_vecs = ws * deltau./ws(:);  % Each row is the "previous point" projected on a different plane. We will have to apply the correct one for the correct points

ptsbm = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^(Nc-1)-1))));  % using binary for the construction of points on a unit square
Lm = deltau^(Nc-1);

ijs = cell2mat(cellfun(@(c) c(:), tausn, 'UniformOutput', false)');  % indices of all points
inds = reshape((1:Nt^Nc)', [repmat(Nt, 1, Nc) ones(Nc==1)]);
S.subs = repmat({':'}, 1, Nc);
S.type = '()';
bpis = cell(Nt,1);
bpjs = cell(Nt,1);
vals = cell(Nt,1);
tic
for ti=1:Nt  % march over the diagonal
    bpis{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,1);
    bpjs{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1));
    vals{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1));
    for ci=1:Nc
        S.subs = repmat({ti:Nt}, 1, Nc);
        S.subs{ci} = ti;
        sinds = unique(reshape(subsref(inds, S), [], 1));
    
        bpis{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1))) = sinds;

        % build local Lagrange interpolant for these points
        cinds = setxor(1:Nc, ci);
        dt_loc = -dt_vecs(ci, cinds);

        n1s_loc = floor(dt_loc/deltau);  % relevant point on the diagonal
        
        xis = (n1s_loc+ptsbm)*deltau;  % points on a cube around previous projection

        Nsfs = abs(prod(xis(end:-1:1,:)-dt_loc,2)'/Lm);  % Lagrangian Shape Functions
        
        ptsb = zeros(2^(Nc-1), Nc);
        ptsb(:, cinds) = n1s_loc+ptsbm+1;  % this +1 is important. Because.

        bpevijs = mod(repmat(ijs(evns(sinds),:), 1, 1, 2^(Nc-1)) + permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the relevant cell
        bpjs{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:) = squeeze(sum((Nt.^(0:Nc-1)).*(bpevijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) 
        vals{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:) = repmat(Nsfs, (Nt-ti+1)^(Nc-1), 1);
    end
    [bpis{ti}, indtis] = unique(bpis{ti});
    bpjs{ti} = bpjs{ti}(indtis,:);
    vals{ti} = vals{ti}(indtis,:);
end
[bpisv, uinds] = unique(cell2mat(bpis));
bpjsv = cell2mat(bpjs);  bpjsv = bpjsv(uinds,:);
valsv = cell2mat(vals);  valsv = valsv(uinds,:);
toc

% Build sparse interpolation matrix
tic
Nmat3 = sparse(repmat(bpisv, 1, 2^(Nc-1)), bpjsv, valsv);
toc
% note to self: Use Nmat as "Dm" is used in scheme 6 in "d_eldrfr1dmar.m"

%% Choose Nmat
switch Nmtype
    case 1
        Nmat = Nmat1;
    case 2
        Nmat = Nmat2;
    case 3
        Nmat = Nmat3;
end

%% Displacement (first harmonic displacement for each component)
hid = eye(Nc);

amps = (1+(1:Nc));
utau = amps(1)*cos(taus{1});
for ti=2:Nc
    utau = utau + amps(ti)*cos(taus{ti});
end
ut = @(t) cos(t(:).*(hfrc*ws(:))')*amps;

%% Elastic Dry Friction Parameters
kt = 1.0;
% muN = sqrt(sum(amps.^2)/2)/2;
switch Nc
    case 1
        muN = 1;
    otherwise
        muN = 2.5;
end

%% Implicit "March"
cstick = @(f) abs(kt*(speye(Nt^Nc)-Nmat)*utau(:) + Nmat*f)<muN;  % stick case
cslip = @(f) abs(kt*(speye(Nt^Nc)-Nmat)*utau(:) + Nmat*f)>=muN;  % slip case

fspr = @(f) f - (kt*(speye(Nt^Nc)-Nmat)*utau(:)+Nmat*f).*cstick(f) - (muN*sign(kt*(speye(Nt^Nc)-Nmat)*utau(:)+Nmat*f)).*cslip(f);
fsp = @(f) deal(f - (kt*(speye(Nt^Nc)-Nmat)*utau(:)+Nmat*f).*cstick(f) - (muN*sign(kt*(speye(Nt^Nc)-Nmat)*utau(:) + Nmat*f)).*cslip(f), ...
    speye(Nt^Nc) - (Nmat).*cstick(f));

tic
opt = struct('Display', true, 'ITMAX', 20);
fsol = NSOLVE(fsp, kt*utau(:), opt);
toc

% figure(Nmtype); surf(reshape(fsol, [repmat(Nt, 1, Nc) ones(1,Nc==1)]))

% Q: DOES FSOLVE HANDLE SPARSE JACOBIANS BETTER THAN A VANILLA IMPLEMENTATION?!

% opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% fsol = fsolve(fsp, kt*utau(:), opt);

% opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display', 'iter');
% fsol = fsolve(fspr, kt*utau(:), opt);

%% Explicit Sequential March with Nmat3
if Nc==2
    Nmat = Nmat3;
    
    fsol_m = zeros(Nt^Nc,1);
    utau_m = utau(:);
    it = 0;
    ress = [];
    while it==0 || ress(end)>eps
        fprev = fsol_m(bpis{1});
        for ti=1:Nt
            ics = bpis{ti};
            
            fsp = kt*(utau_m(ics) - Nmat(ics, :)*utau_m) + Nmat(ics, :)*fsol_m; % stick prediction
    
            fsol_m(ics) = fsp.*(abs(fsp)<muN) + muN*sign(fsp).*(abs(fsp)>muN);
    
    %         clf()
    %         plot(taus{1}, taus{2}, 'k.'); hold on
    %         plot(taus{1}(ics), taus{2}(ics), 'ro')
    %         keyboard
        end
        it = it+1;
        ress = [ress; vecnorm(fprev-fsol_m(bpis{1}))];
        fprintf('%d %e\n', it, ress(it))
    
        fmx = reshape(fsol_m, [repmat(Nt, 1, Nc) ones(1,Nc==1)]);  % f in matrix form (in tau space)
        
        repinds = mat2cell(repmat([1:Nt 1]', 1, Nc), Nt+1, ones(1, Nc));
        fmxp = fmx(repinds{:});
        
        figure(3)
        clf()
        set(gcf, 'Color', 'white')
        surf(taus{1}, taus{2}, fmx);
        colormap(jet)
        tics = (0:.25:1);
        ticls = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
        set(gca, 'XTick', tics*2*pi)
        set(gca, 'XTickLabel', ticls)
        set(gca, 'YTick', tics*2*pi)
        set(gca, 'YTickLabel', ticls)
        xlabel('Time $\tau_1$', 'Rotation', 28)
        ylabel('Time $\tau_2$', 'Rotation', -38)
        zlabel('Force (N)')
        if plotfigs
            export_fig(sprintf('./FIGS/F_marchsurf_%d.png', it), '-dpng')
        end

    %     plot(fmx, '.-')
        figure(200)
        clf()
        semilogy(ress, '.-')
        grid on
    
        keyboard
    end
end
figure(200)
clf()
set(gcf, 'Color', 'white')
semilogy(ress, 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'w')
grid on
xlabel('Iterations')
ylabel('First front deviation ($L_2$)')
if plotfigs
    export_fig('./FIGS/F_marchconv.png', '-dpng')
end

% %%
% figure(100)
% clf()
% plot(taus{1}, taus{2}, 'k.'); hold on
% plot(taus{1}(bpis{3}), taus{2}(bpis{3}), 'ro')
% plot(Nmat3(bpis{3},:)*taus{1}(:), Nmat3(bpis{3},:)*taus{2}(:), 'b*')
% 
% pti = 3;
% plot(taus{1}(bpis{3}(pti)), taus{2}(bpis{3}(pti)), 'ro', 'MarkerFaceColor', 'r')
% plot(Nmat3(bpis{3}(pti),:)*taus{1}(:), Nmat3(bpis{3}(pti),:)*taus{2}(:), 'bo', 'MarkerFaceColor', 'b')
% 
% plot(taus{1}(bpis{3})'+ws(1)*[-1 1]'*10, taus{2}(bpis{3})'+ws(2)*[-1 1]'*10, 'g-')
% 
% xlim([-deltau 2*pi])
% ylim([-deltau 2*pi])
% 
% axis equal
% axis off

%% Matrix Sparsity
figure(4)
clf()
spy(Nmat3, 30/Nc)
set(gca, 'XTick', fix(linspace(1, Nt^Nc, Nt)))
set(gca, 'YTick', fix(linspace(1, Nt^Nc, Nt)))
xlim([1 Nt^Nc]); ylim([1 Nt^Nc])
grid on
set(gcf, 'Color', 'white')
if plotfigs
    export_fig(sprintf('./FIGS/F_Nmat3_%d.png', Nc), '-dpng')
end