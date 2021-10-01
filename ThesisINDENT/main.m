clc
clear all
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

segi  = 4;
reso = 3.7019;  % Micrometers

Aspstruc = struct('center', [1 2 3]', 'rads', [1 2 3]', ...
    'csys', eye(3), 'Bmat', eye(3), 'Lvec', [1 2 3]', 'Csca', 0);

%% Read Data
dat = dlmread(sprintf('./DATS/R3AHighMag_Subset_%d_Nds_81_by_81.dat', segi))*1e6;  % Microns
Nq = size(dat, 1);
[xq, yq] = meshgrid((0:Nq-1)*reso, (0:Nq-1)*reso);

% Nq2 = 256;
% [xq2, yq2] = meshgrid((0:Nq2-1)*reso*Nq/Nq2, (0:Nq2-1)*reso*Nq/Nq2);
% dat2 = interp2(xq, yq, dat, xq2, yq2);

%% Image Processing
se = strel('square', 2);
dat_o = imerode(dat, se);
dat_o = imdilate(dat_o, se);

[zqe, ~, Gv, Gh] = edge(dat_o);

Gradz = sqrt(Gv.^2+Gh.^2);

dat_s = sort(dat_o(:), 'descend');
Sz = mean(dat_s(1:5)-dat_s(end:-1:end-4));  % Sz(ISO): 10-point height

% WG = watershed(Gradz.*(Gradz>0.005*Sz));
WG = watershed(Gradz);
%% Detect Asperities
count = 0;

% figure(1)
% clf()
% imagesc(WG)
Asperities = repmat(Aspstruc, max(WG(:)), 1);
for i=1:max(WG(:))
    [i1, j1] = find(WG==i);
    ij1 = [i1 j1];
    % Boundaries
    ijbds1 = unique([ij1-[1 0]; ij1-[0 1]; ij1+[1 0]; ij1+[0 1];], 'rows');
    ijbds1 = ijbds1(ijbds1(:,1)<=Nq & ijbds1(:,1)>0 & ijbds1(:,2)<=Nq & ijbds1(:,2)>0, :);
    ijbds1 = setdiff(ijbds1, ij1, 'rows');
    % Indices
    ind1 = sub2ind([Nq Nq], ij1(:,1), ij1(:,2));
    bdind1 = sub2ind([Nq Nq], ijbds1(:,1), ijbds1(:,2));
    
    isasp = (((mean(dat_o(ind1))-mean(dat_o(bdind1)))/abs(mean(dat_o(ind1))))>0.10 & ...
        length([ind1; bdind1])>=9);% & ...
        %~any(i1==1 | i1==Nq | j1==1 | j1==Nq));
%         if mean(i1)<10 && mean(j1)>50 && mean(j1)<60
%             keyboard
%         end
        
    % Check for asperity
    if isasp %% && sum(all(dat_o(ind1) > dat_o(bdind1)', 2))>length(ind1)/8       
       % Fit Conic Section
       bdind2 = bdind1(dat(bdind1)<min(dat(ind1)));
       inds = [ind1; bdind1];
       xyz = [xq(inds) yq(inds) dat(inds)];
       try
           % Generic fit (hyperboloid, paraboloid, ellipsoid)
%             [center, axesq, csys, Serr, Bmat, Lvec, Csca] = ELLIPSOID3D(xyz,1e-10);

           % Restrict to ellipsoid aligned on z (x-y orthotropy allowed)
           [center, axesq, csys, Serr, Bmat, Lvec, Csca] = ELLIPSOID3D_zfix(xyz,1e-10);
       catch me  % Can't avoid single sheet
%            disp(me)
           continue;
       end    

        %% Asperity Model
%         hold on; plot(j1, i1, 'ko', 'MarkerFaceColor', 'k'); 
        
        count = count+1;
        Asperities(count).center = center;
        Asperities(count).rads   = sqrt(axesq);
        Asperities(count).csys   = csys;
        Asperities(count).Bmat   = Bmat;
        Asperities(count).Lvec   = Lvec;
        Asperities(count).Csca   = Csca;
    end
end
Asperities = Asperities(1:count);

%% Identification
cs = [Asperities.center];
as = [Asperities.rads];
zs = cs(3,:)+as(3,:);

[pvals, zvals] = ecdf(zs);
[lx0l, stats] = robustfit(-zvals(1:end-1), log(1-pvals(1:end-1)));

% figure(2)
% clf()
% plot(zvals, pvals); hold on
% plot(zvals, 1-exp(-lx0l(2)*(zvals-lx0l(1)/lx0l(2))), 'r-.')

%% Contact Model
Nasps = count;
lam = lx0l(2);
Rad = mean([Asperities.rads]');  Rad = Rad(3);

E = 200000;
nu = 0.3;

Estar = E/(1-nu^2);

cn = Estar*sqrt(pi*Rad/lam^3)*Nasps/4;  %% This scaling is adjusted since the estimate of the number of asperities is not relied upon, so we can slide the curve on the horizontal axis

%% Plots
figure(1)
clf()
set(gcf, 'Color', 'white')

aa = gobjects(5, 1);

k = 1;
fudat = dlmread(sprintf('./DATS/ref.dat'));
aa(k) = plot(fudat(:,2), -fudat(:,1)*1e-6, '-', 'LineWidth', 2); hold on
legend(aa(k), 'Flat Reference')
k = k+1;
for mk=[21 41 81]  % Different meshes
    fudat = dlmread(sprintf('./DATS/S4_M%d.dat', mk));
    aa(k) = plot(fudat(:,2), -fudat(:,1)*1e-6, '-', 'LineWidth', 2); hold on
    
    legend(aa(k), sprintf('EP-FE, Mesh: $%d \\times %d$', mk, mk), 'interpreter', 'latex')
    k = k+1;
end
ds = linspace(-10, 25, 1000);
aa(k) = plot(ds, cn.*exp(lam*(ds))*1e-6, 'k-.', 'LineWidth', 3);
legend(aa(k), 'Rough Contact Model')
ylim([0 1.15e2])
xlim([-10 20])

ll = legend(aa(1:end), 'Location', 'southeast');
xlabel('Planar Approach ($\mu\,m$)')
ylabel('Normal Force (N)')

export_fig('./FIGS/FE_RC_results.eps', '-depsc')