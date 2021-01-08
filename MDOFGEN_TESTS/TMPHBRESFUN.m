function [R, dRdU, dRdw] = TMPHBRESFUN(Uw, M, C, K, Fl, kt, kn, mu, Nt, h)
    U = Uw(1:end-1);
    w = Uw(end);
    
    Nhc = sum((h==0)+2*(h~=0));
    
    [E, dEdw] = HARMONICSTIFFNESS(M, C, K, w, h);
    
    % Manual friction model evaluation
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);
    
    uxt = TIMESERIES_DERIV(Nt, h, U(1:2:end), 0);
    unt = TIMESERIES_DERIV(Nt, h, U(2:2:end), 0);
    
    fxt = zeros(Nt, 1);
    fnt = zeros(Nt, 1);
    
    jxx = zeros(Nt, Nhc);
    jxn = zeros(Nt, Nhc);    
    jnn = zeros(Nt, Nhc);
    
    % NORMAL CONTACT
    fnt = max(kn*unt, 0);
    jnn = kn.*(fnt~=0).*cst;
    
    % TANGENTIAL CONTACT
% 	fxt = kt*uxt;
%     jxx = kt*cst;
    while (1)
        fprev = fxt(end);
        for ti=1:Nt
            tm1 = mod(ti-1-1, Nt)+1;

            if fnt(ti)~=0  % in contact
                fslip = abs(mu*fnt(ti));  % slip force
                fsp = kt*(uxt(ti)-uxt(tm1)) + fxt(tm1);  % stick prediction

                if abs(fsp)<fslip  % stick
                    fxt(ti) = fsp;

                    jxx(ti, :) = kt*(cst(ti,:)-cst(tm1,:)) + jxx(tm1, :);
                    jxn(ti, :) = jxn(tm1, :);
                else  % slip
                    fxt(ti) = fslip*sign(fsp);

                    jxx(ti, :) = 0;
                    jxn(ti, :) = mu*kn*sign(fsp)*cst(ti, :);

    %                 jxn(ti, :) = mu*kn*cst(ti, :);
                end
                else  % separated
                fxt(ti) = 0;
                jxx(ti, :) = 0;
                jxn(ti, :) = 0;
            end
        end
        if abs(fprev-fxt(end))<1e-6
            break;
        end
%         disp('reping')
    end
    
%     fxt = kt*uxt;
%     fnt = kn*unt;
%     jxx = kt*cst;
%     jnn = kn*cst;
%     jxn = zeros(Nhc);
    
    FN = GETFOURIERCOEFF(h, fnt);
    FX = GETFOURIERCOEFF(h, fxt);
    JXX = GETFOURIERCOEFF(h, jxx);
    JXN = GETFOURIERCOEFF(h, jxn);
    JNN = GETFOURIERCOEFF(h, jnn);
    % Manual friction model evaluation
    
    % Residue
    Fnl = zeros(Nhc*2,1);
    Fnl(1:2:end) = FX;
    Fnl(2:2:end) = FN;
    
    Jnl = zeros(Nhc*2);
    Jnl(1:2:end, 1:2:end) = JXX;
    Jnl(1:2:end, 2:2:end) = JXN;
	Jnl(2:2:end, 2:2:end) = JNN;
    
	R = E*U + Fnl - Fl;
    dRdU = E + Jnl;
    dRdw = dEdw*U;
    dRdUw = [dRdU dRdw];
end