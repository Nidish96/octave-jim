function [R, dRdU, FNL] = QPHBRESFUN(m, Up, ws, Fl, h, Nt, tol, varargin)
%QPHBRESFUN 
%
%   USAGE: 
%       [R, dRdU, dRdw, FNL] = QPHBRESFUN(m, Up, ws, Fl, h, Nt, tol)
%   INPUTS:
%       MDOFGEN class
%       Up
%       Fl
%       h
%       Nt
%       tol 
%   OUTPUTS:
%       R
%       dRdU
%       dRdw
%       FNL
    
    Nc = size(h,2);  % Number of components
    Nh = size(h,1);  % Number of harmonics
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));  % Number of harmonic coefficients
    
    Ws = Up(end)*ws;  % real frequency is ws scaled by Up(end)

    E = QPHARMONICSTIFFNESS(m.M, m.C, m.K, Ws, h);  % Harmonic Stiffness
    
    tau = linspace(0, 2*pi, Nt+1); tau(end) = [];
    taus = cell(Nc, 1);
    [taus{:}] = ndgrid(tau);
    t = zeros(repmat(Nt, 1, Nc));
    for ic=1:Nc
        t = t+taus{ic}*Ws(ic);
    end
    
    D1 = QPHARMONICSTIFFNESS(0, 1, 0, Ws, h);  % time derivative matrix
    
    cst = QPTIMETRANS(eye(Nhc), h, Nt);  % basis functions
    sct = QPTIMETRANS(D1, h, Nt);  % basis function derivatives
    
    FNL = zeros(m.Ndofs*Nhc, 1);
    dFNL = zeros(m.Ndofs*Nhc);
    for ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(Up(1:end-1), m.Ndofs, Nhc))';  % Nhc x Ndnl
        Ndnl = size(m.NLTs(ni).L, 1);

        unlt   = QPTIMETRANS(Unl, h, Nt);
        unldot = QPTIMETRANS(D1*Unl, h, Nt);
        
        % INSTANTANEOUS FORCING
        [ft, dfdu, dfdud] = m.NLTs(ni).func(taus, unlt, unldot);
        
        F = QPFOURIERCOEFF(ft, h);
        J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
        dFdU = reshape(QPFOURIERCOEFF(reshape(dfdu.*permute(cst, [1, 3, 2]), Nt^Nc, Ndnl*Nhc), h), ...
            Nhc, Ndnl, Nhc);
        
        for di=1:Ndnl
            J(di:Ndnl:end, di:Ndnl:end) = dFdU(:, di, :);
        end
        
        if m.NLTs(ni).type<=5  % Self adjoint forcing
            FNL  = FNL + reshape(m.NLTs(ni).L'*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).L')*J*kron(eye(Nhc), m.NLTs(ni).L);
        else  % Non-self adjoint forcing
            FNL  = FNL + reshape(m.NLTs(ni).Lf*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).Lf)*J*kron(eye(Nhc), m.NLTs(ni).L);                
        end
    end
    
    % Residue
    if length(m.Rsc)~=length(Fl)
        m.Rsc = (1/max(abs(Fl)))*ones(length(Fl),1);
    end
    R = [E*Up(1:end-1) + FNL - Fl].*m.Rsc;
    dRdU = (E+dFNL).*m.Rsc;
%     dRdw = (dEdw*Up(1:end-1)).*m.Rsc;

    % All Gradient terms in one matrix
%     if ~isempty(varargin)
%         dRdU = [dRdU dRdw];
%     end
end