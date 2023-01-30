function [R, dRdUwxi, dRdlA, dRdash, FNL] = EQPMCRESFUN(m, UwxiA, ash, Fls, h, Nt, tol, varargin)
%QPHBRESFUN 
%   NOTE: h MUST contain (0,1) & (1,0) components
%
%       VELOCITY DEPENDENT NONLINEARITIES NOT IMPLEMENTED YET!
%
%   USAGE: 
%       [R, dRdU, dRdw, FNL] = EQPMCRESFUN(m, UwxiA, ash, Fls, h, Nt, tol)
%   INPUTS:
%       MDOFGEN class
%       UwxiA       : (Nhc*Nd+Nc+Nc+1, 1)
%       ash         : (Nc,1) Amplitude vector 
%       Fls         : (Nhc*Nd, Nc) Phase constraint matrix
%       h           : (Nh, Nc)
%       Nt          : (1,1)
%       tol         : (1,1)
%   OUTPUTS:
%       R           : (Nhc*Nd+2*Nc, 1)
%       dRdUwxi     : (Nhc*Nd+2*Nc, Nhc*Nd+2*Nc)
%       dRdlA       : (Nhc*Nd+2*Nc, 1)
%       dRdlAs      : (Nhc*Nd+2*Nc, Nc)
%       FNL         : (Nhc*Nd, 1)

    % Uwxia: [Uh; ws; xis; log(A)]

    Nc = size(h,2);  % Number of components
    Nh = size(h,1);  % Number of harmonics
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));  % Number of harmonic coefficients
    
    xis = UwxiA(end-Nc:end-1);
    ws = UwxiA(end-2*Nc:end-Nc-1);
    
    lA = UwxiA(end);
    A = 10^lA;
    dAdlA = A*log(10);
    
    h0 = double(all(h(1,:)==0));
    Asc = kron([ones(h0,1); A*ones(Nhc-h0,1)], ones(m.Ndofs,1));
    dAscdA = kron([zeros(h0,1); ones(Nhc-h0,1)], ones(m.Ndofs,1));
    
    As = A*ash;  % List of modal amplitudes
    
    Uh = Asc.*UwxiA(1:end-2*Nc-1);  % Harmonic Coefficients of Mode Shape

    Cg = m.CAUGHYMATS(Nc, 0);  % Caughy Matrices
    xiM = sum(reshape(xis, [1 1 Nc]).*Cg, 3);
    [E, dEdws] = QPHARMONICSTIFFNESS(m.M, m.C-xiM, m.K, ws, h);  % Harmonic Stiffness
    dEdxi = zeros(m.Ndofs*Nhc, m.Ndofs*Nhc, Nc);  % Maybe using cells of sparse matrices is better for storage
    for ci=1:Nc  % Can be replaced with pagefun
        dEdxi(:, :, ci) = QPHARMONICSTIFFNESS(m.M*0, -Cg(:, :, ci), m.M*0, ws, h);
    end
    
    [FNL, dFNL, dFNLdw] = m.QPNLEVAL([Uh; ws(:)], h, Nt, tol);
    
    % Residue        
    acons = zeros(Nc, 1);
    d_acons = zeros(Nc, Nhc*m.Ndofs+2*Nc);
    d_acons_lA = zeros(Nc, 1);
    d_acons_lAs = zeros(Nc, Nc);  % A=As^2/ash^2
    d_acons_ash = zeros(Nc, Nc);
    for ci=1:Nc
        inds = setxor(1:Nc, ci);
        
        hid = find(all(h(:, [ci inds])==[1 zeros(1, Nc-1)], 2));
        ampis = m.Ndofs*h0+(hid-1-h0)*2*m.Ndofs + (1:m.Ndofs);  % 1 x Nd
        
        acons(ci) = (Uh(ampis)'*m.M*Uh(ampis) + Uh(m.Ndofs+ampis)'*m.M*Uh(m.Ndofs+ampis)) - As(ci)^2;
        d_acons(ci, ampis) = 2*Uh(ampis)'*m.M*A;
        d_acons(ci, m.Ndofs+ampis) = 2*Uh(m.Ndofs+ampis)'*m.M*A;
        d_acons_lA(ci) = (2*(Uh(ampis)'*m.M*Uh(ampis) + Uh(m.Ndofs+ampis)'*m.M*Uh(m.Ndofs+ampis))/A - 2*As(ci)*ash(ci)^2)*dAdlA;
    
        d_acons_ash(ci, ci) = 2*A^2*ash(ci);
    end
    
    R = [E*Uh+FNL;    % balance equations     (Nhc*Nd)
        acons;          % amplitude constrains  (Nc)
        Fls'*(Uh./Asc)];      % phase constraints     (Nc)
    dRwx = zeros(Nhc*m.Ndofs, 2*Nc);
    for ci=1:Nc
        dRwx(:, ci) = dEdws(:, :, ci)*Uh;  % Needs an additional dFNLdws term for velocity dependent nls
        dRwx(:, Nc+ci) = dEdxi(:, :, ci)*Uh;
    end
    dRwx(:, 1:Nc) = dRwx(:, 1:Nc) + dFNLdw;
    
    dRdUwxi = [(E+dFNL)*diag(Asc), dRwx;
        d_acons;
        Fls', zeros(Nc, 2*Nc)];
    
    dRdlA = [(E+dFNL)*(Uh./Asc).*dAscdA*dAdlA;
        d_acons_lA;
        zeros(Nc, 1)];
    
    dRdash = [zeros(Nhc*m.Ndofs, 2);
              d_acons_ash;
              zeros(Nc, 2)];

    % All Gradient terms in one matrix
    if ~isempty(varargin)
        dRdUwxi = [dRdUwxi dRdlA];
    end
end

%%
function [res, dresdf, dresdu] = QPMARCHRESFUN(ft, unlt, func, Nmat)
%Residue function for 
    [fnl, ~, dfnldfp, dfnlduc, dfnldup] = func(unlt, Nmat*unlt, Nmat*ft);
    
    res = ft-fnl;
    dresdf = sparse(eye(length(ft))-diag(dfnldfp)*Nmat);
    dresdu = sparse(-diag(dfnlduc) - diag(dfnldup)*Nmat);
%             [ft(ti,:), dfdu(ti,:,:,:)] = ...
%                 m.NLTs(ni).func(t(ti), unlt(ti,:), h, t(tm1), ...
%                 unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:,:));    
end
%%
function Nmat = CONSTRUCTNMAT(ws, Nc, Nt, varargin)  % Need to think about w-gradients
    if nargin==3
        Nmtype = 1;
    else
        Nmtype = varargin{1};
    end
    taus = cell(Nc, 1);
    [taus{:}] = ndgrid(1:Nt);    
    switch Nmtype
        case 1
            deltau = 2*pi/Nt; 
            dt_vec = ws*deltau/vecnorm(ws);  % vector corresponding to deltatau amplitude in real time dxn

            ptsb = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^Nc-1))));  % using binary for the construction of points on a unit square
            oppi = (2^Nc):-1:1;    % diagonally opposite points are retrieved using the binary inverses
            xis = ptsb*(-deltau);  % coordinates in "tau" space relative to origin

            Lm = deltau^Nc;                                 % Lebesgue Measure of each cell in tau space
            Nsf = prod(abs(xis(oppi,:)-(-dt_vec)), 2)'/Lm;  % Lagrange Shape Functions to interpolate previous point in the cell diagonally behind (in tau space)

            ijs = fix(cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)'));  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, 2^Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;   % indices of points forming the cell that is diagonally behind
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);         % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

            % Build sparse interpolation matrix
            Nmat = sparse(repmat((1:Nt^Nc)', 1, 2^Nc), evns, repmat(Nsf, Nt^Nc, 1));
        case 2
            ptsb = eye(Nc);

            Nsf = ws/sum(ws);  % Shape functions constructed using omegas

            ijs = cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)');  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the FD stencil
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points from the stencil

            % Build sparse interpolation matrix
            Nmat = sparse(repmat((1:Nt^Nc)', 1, Nc), evns, repmat(Nsf, Nt^Nc, 1));
    end
end
