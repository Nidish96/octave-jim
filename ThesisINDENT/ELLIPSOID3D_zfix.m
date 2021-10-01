function [center, axes, csys, Serr, Bmat, Lvec, Csca] = ELLIPSOID3D_zfix(xyz, tol)
%ELLIPSOID3D_zfix fits 3D data to 3D ellipsoid (rotated) where only
%rotations about the z axis are permitted
%
%   USAGE:
%       [center, axes, csys] = ELLIPSOID3D_zfix(xyz)
%   INPUTS:
%       xyz     : (Npt, 3)
%   OUTPUTS:
%       center  : (3, 1) center
%       axes    : (3, 1) semi-minor axes lengths
%       csys    : (3, 3) coordinate transform matrix (affine)

    if size(xyz,1)==3 && size(xyz,2)~=3
        xyz = xyz';
    end
    
    if nargin==1
        tol = eps;
    end
    
	% Scale Data
    mxyz = mean(xyz);
%     mxyz = [1 1 1];
    
    Tm = diag(mxyz);
    Tmi = diag(1./mxyz);
    xyz = xyz*Tmi;
    
    [~, S, V] = svd([xyz.^2 xyz(:,1).*xyz(:,2) xyz ones(size(xyz,1),1)]);
    S = diag(S);
    Serr = S(end)/S(1);
    acofs = V(:, end);

    % Assemble Matrices (rescaled)
    Bmat = Tmi'*[acofs(1) acofs(4)/2 0;
        acofs(4)/2 acofs(2) 0;
        0 0 acofs(3)]*Tmi;
    Lvec = Tmi*acofs(5:7);
    Csca = acofs(8);

    center = -Bmat\Lvec/2;

    [bU, bS, bV] = svd(Bmat);
    sc = center'*Bmat*center-Csca;
    eV = bU;
    eD = diag(bU'*Bmat*bU)/sc;
    
    si = [1 2 3]; 
    [~, si(1)] = max(abs(eV(1,:)));
    sim = setdiff([1 2 3], si(1));
    [~, si(2)] = max(abs(eV(2,sim))); si(2) = sim(si(2));
    si(3) = setdiff([1 2 3], si(1:2));

    eD = eD(si); eV = eV(:, si);
    
    eV = eV*diag(sign(diag(eV)));
    axes = 1./eD;  % axes^2 of ellipse (or hyperboloid/praboloid if negative)
    axes(~isfinite(axes)) = 0.0;
    csys = eV;
    
    if sum(axes>0)~=3
        error('Not an ellipse')
    end
    
    % Matrices for return
    Bmat = Bmat/sc;
    Lvec = Lvec/sc;
    Csca = Csca/sc;
end