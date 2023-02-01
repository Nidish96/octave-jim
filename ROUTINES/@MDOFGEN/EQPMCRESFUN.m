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
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));  % Number of harmonic coefficients
    
    xis = UwxiA(end-Nc:end-1);
    ws = UwxiA(end-2*Nc:end-Nc-1);
    
    lA = UwxiA(end);
    A = 10^lA;
    dAdlA = A*log(10);
    
    [rinds0,zinds,hinds,rinds,iinds] = HINDS(m.Ndofs, h);
    Asc = ones(Nhc*m.Ndofs,1);  Asc([rinds; iinds]) = A;
    dAscdA = zeros(Nhc*m.Ndofs,1); dAscdA([rinds; iinds]) = 1;

    h0 = double(all(h(1,:)==0));

    As = A*ash;  % List of modal amplitudes
    
    Uh = Asc.*UwxiA(1:end-2*Nc-1);  % Harmonic Coefficients of Mode Shape

    Cg = m.CAUGHYMATS(Nc, 0);  % Caughy Matrices
    xiM = sum(reshape(xis, [1 1 Nc]).*Cg, 3);
    [E, dEdws] = QPHARMONICSTIFFNESS(m.M, m.C-xiM, m.K, ws, h);  % Harmonic Stiffness
    dEdxi = zeros(m.Ndofs*Nhc, m.Ndofs*Nhc, Nc);  % Maybe using cells of sparse matrices is better for storage
    for ci=1:Nc  % Can be replaced with pagefun
        dEdxi(:, :, ci) = QPHARMONICSTIFFNESS(zeros(m.Ndofs), -Cg(:, :, ci), zeros(m.Ndofs), ws, h);
    end
    if any(isnan(ws))
        keyboard
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
        dRwx(:, ci) = dEdws(:, :, ci)*Uh + dFNLdw(:,ci);  % Needs an additional dFNLdws term for velocity dependent nls
        dRwx(:, Nc+ci) = dEdxi(:, :, ci)*Uh;
    end
    
    dRdUwxi = [(E+dFNL)*diag(Asc), dRwx;
        d_acons;
        Fls', zeros(Nc, 2*Nc)];
    
    dRdlA = [(E+dFNL)*(UwxiA(1:end-2*Nc-1).*dAscdA)*dAdlA;
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