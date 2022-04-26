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
        
        if mod(m.NLTs(ni).type, 2)==0  % INSTANTANEOUS FORCING
            [ft, dfdu, dfdud] = m.NLTs(ni).func(taus, unlt, unldot);
        
            F = QPFOURIERCOEFF(ft, h);
            J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
            dFdU = reshape(QPFOURIERCOEFF(reshape(dfdu.*permute(cst, [1, 3, 2]), Nt^Nc, Ndnl*Nhc), h) + ...
                QPFOURIERCOEFF(reshape(dfdud.*permute(sct, [1, 3, 2]), Nt^Nc, Ndnl*Nhc), h), ...
                Nhc, Ndnl, Nhc);
        
            for di=1:Ndnl
                J(di:Ndnl:end, di:Ndnl:end) = dFdU(:, di, :);
            end
        else % HYSTERETIC FORCING
            Nnlf = size(m.NLTs(ni).L, 1);
            if ~isempty(m.NLTs(ni).Lf)
                Nnlf = size(m.NLTs(ni).Lf, 2);
            end
            if Nnlf>1 || size(unlt, 2)>1
                error('Not implemented for multi-DOF dependent multi-force nonlinearities yet')
            end

            % Construct Nmat
            Nmat = CONSTRUCTNMAT(Ws, Nc, Nt, m.NLTs(ni).qptype);
            ft = ones(Nt^Nc, Nnlf);
            
            opt = struct('Display', false);
            ft = NSOLVE(@(ft) QPMARCHRESFUN(ft, unlt, ...
                @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
                Nmat), ft, opt);
%             opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'off');
%             ft = fsolve(@(ft) QPMARCHRESFUN(ft, unlt, ...
%                 @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
%                 Nmat), ft, opt);

            [~, dresdf, dresdu] = QPMARCHRESFUN(ft, unlt, ...
                @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
                Nmat);
            dfdut = -dresdf\dresdu;
            
            F = QPFOURIERCOEFF(ft, h);
            J = QPFOURIERCOEFF(dfdut*cst, h);
%             J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
%             for di=1:Ndnl
%                 for dj=1:Ndnl
%                     tmp = squeeze(dfdu(:, di, dj, :));
%                     if ~isempty(find(tmp~=0, 1))
%                         J(di:Ndnl:end, dj:Ndnl:end) = ...
%                             GETFOURIERCOEFF(h, tmp);
%                     end
%                 end
%             end
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