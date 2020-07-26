function [R, dRdU, dRdw, FNL] = HBRESFUN(m, Uw, Fl, h, Nt, tol)

    Nhc = sum((h==0)+2*(h~=0));
    
    w = Uw(end);
    
    [E, dEdw] = HARMONICSTIFFNESS(m.M, m.C, m.K, w, h);
    
    t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  
    sct = w*TIMESERIES_DERIV(Nt, h, eye(Nhc), 1);
    
    FNL = zeros(m.Ndofs*Nhc, 1);
    dFNL = zeros(m.Ndofs*Nhc);
    for ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(Uw(1:end-1), m.Ndofs, Nhc))';  % Nhc x Ndnl
        Ndnl = size(m.NLTs(ni).L, 1);
        
        unlt = TIMESERIES_DERIV(Nt, h, Unl, 0);  % Nt x Ndnl
        unldot = w*TIMESERIES_DERIV(Nt, h, Unl, 1);  % Nt x Ndnl
        
        if mod(m.NLTs(ni).type-1, 3)==0  % Instantaneous force
            [ft, dfdu, dfdud] = m.NLTs(ni).func(t, unlt, unldot);
            % (Nt,Ndnl); (Nt,Ndnl); (Nt,Ndnl) (point-wise)
        else  % Hysteretic force
            ft = zeros(Nt, Nd);
            dfdu = zeros(Nt, Nd, Nhc);
            
            fprev = ft(end, :);
            its = 0;
            
            while abs(fprev-ft(end, :))>tol || its==0
                for ti=1:Nt
                    tm1 = mod(ti-2, Nt)+1;
                    [ft(ti,:), dfdu(ti,:,:)] = ...
                        m.NLTs(ni).func(t(ti), unlt(ti,:), t(mt1), ...
                        unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:), h);
                end
            end
        end
        F = GETFOURIERCOEFF(h, ft);
        dFdU = reshape(GETFOURIERCOEFF(h, reshape(dfdu.*permute(cst, [1, 3, 2]), Nt, Ndnl*Nhc)),...
            Nhc, Ndnl, Nhc);
        J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
        for di=1:Ndnl
            J(di:size(m.NLTs(ni).L,1):end, di:size(m.NLTs(ni).L,1):end) = dFdU(:, di, :);
        end
        
        if m.NLTs(ni).type<=5  % Self adjoint forcing
            FNL = FNL + reshape(m.NLTs(ni).L'*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).L')*J*kron(eye(Nhc), m.NLTs(ni).L);
	    else  % Non-self adjoint forcing
            FNL = FNL + reshape(m.NLTs(ni).Lf*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).Lf)*J*kron(eye(Nhc), m.NLTs(ni).L);                
        end
    end
    
    % Residue
    R = [E*Uw(1:end-1) + FNL - Fl];
    dRdU = E+dFNL;
    dRdw = dEdw*Uw(1:end-1);
end