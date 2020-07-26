function [R, dRdUwx, dRda] = EPMCRESFUN(m, Uwxa, Fl, h, Nt, tol)

    Nhc = sum((h==0)+2*(h~=0));

    la = Uwxa(end);
    A = 10^la;
    dAdla = A*log(10);
    
    xi = Uwxa(end-1);
    w = Uwxa(end-2);
    
    [E, dEdw] = HARMONICSTIFFNESS(m.M, m.C-xi*m.M, m.K, w, h);
    dEdxi = HARMONICSTIFFNESS(m.M*0, -m.M, m.M*0, w, h);
    
    t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  
    sct = w*TIMESERIES_DERIV(Nt, h, eye(Nhc), 1);
    
    FNL = zeros(m.Ndofs*Nhc, 1);
    dFNL = zeros(m.Ndofs*Nhc);
    for ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(A*Uwxa(1:end-3), m.Ndofs, Nhc))';  % Nhc x Ndnl
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
    Nd = m.Ndofs;    
%     R = [E*(A*Uwxa(1:end-3)) + FNL;
%         Uwxa(Nd*(h(1)==0)+(1:Nd))'*m.M*Uwxa(Nd*(h(1)==0)+(1:Nd))+Uwxa(Nd*(h(1)==0)+Nd+(1:Nd))'*m.M*Uwxa(Nd*(h(1)==0)+Nd+(1:Nd))-1;
%         Fl'*Uwxa(1:end-3)];
%     dRdUwx = [(E+dFNL)*A, dEdw*(A*Uwxa(1:end-3)), dEdxi*(A*Uwxa(1:end-3));
%         zeros(1, Nd*(h(1)==0)), Uwxa(Nd*(h(1)==0)+(1:2*Nd))'*blkdiag(M,M), zeros(1, Nd*(Nhc-2-(h(1)==0))), 0, 0;
%         Fl', 0, 0];
%     dRda = [(E+dFNL)*Uwxa(1:end-3)*dAdla;
%         0;
%         0];
    
    R = [E*(A*Uwxa(1:end-3))+FNL;
        Uwxa(1:end-3)'*kron(diag([0 ones(1,Nhc-1)]), m.M)*Uwxa(1:end-3)-1;
        Fl'*Uwxa(1:end-3)];
    dRdUwx = [(E+dFNL)*A, dEdw*(A*Uwxa(1:end-3)), dEdxi*(A*Uwxa(1:end-3));
        2*Uwxa(1:end-3)'*kron(diag([0 ones(1,Nhc-1)]), m.M), 0, 0;
        Fl', 0, 0];
    dRda = [(E+dFNL)*Uwxa(1:end-3)*dAdla;
        0; 
        0];
end