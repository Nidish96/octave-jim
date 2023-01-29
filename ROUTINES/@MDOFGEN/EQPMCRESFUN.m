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
    Ws = UwxiA(end-2*Nc:end-Nc-1);
    
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
    [E, dEdws] = QPHARMONICSTIFFNESS(m.M, m.C-xiM, m.K, Ws, h);  % Harmonic Stiffness
    dEdxi = zeros(m.Ndofs*Nhc, m.Ndofs*Nhc, Nc);  % Maybe using cells of sparse matrices is better for storage
    for ci=1:Nc  % Can be replaced with pagefun
        dEdxi(:, :, ci) = QPHARMONICSTIFFNESS(m.M*0, -Cg(:, :, ci), m.M*0, Ws, h);
    end
    
    tau = linspace(0, 2*pi, Nt+1); tau(end) = [];
    taus = cell(Nc, 1);
    [taus{:}] = ndgrid(tau);
    t = zeros(repmat(Nt, 1, Nc));
    for ic=1:Nc
        t = t+taus{ic}*Ws(ic);
    end
    
    [D1, dD1dws] = QPHARMONICSTIFFNESS(0, 1, 0, Ws, h);  % time derivative matrix
    
    cst = QPTIMETRANS(eye(Nhc), h, Nt);  % basis functions
    sct = QPTIMETRANS(D1, h, Nt);  % basis function derivatives
    
    dsct_ws = QPTIMETRANS(reshape(dD1dws, Nhc, Nhc*Nc), h, Nt);  % Nt^Nc x Nhc*Nc ("pages" next to each other)
    
    FNL = zeros(m.Ndofs*Nhc, 1);
    dFNL = zeros(m.Ndofs*Nhc);
    for ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(Uh, m.Ndofs, Nhc))';  % Nhc x Ndnl
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
%             dFdws = reshape(QPFOURIERCOEFF(reshape(dfdud.*permute(dsct_ws, [1, 3, 2]), Nt^Nc, Ndnl*Nhc*Nc), h), ...
%                 Nhc, Ndnl, Nhc, Nc);
            
%             Jw = zeros([size(J) Nc]);
            for di=1:Ndnl
                J(di:Ndnl:end, di:Ndnl:end) = dFdU(:, di, :);
%                 Jw(di:Ndnl:end, di:Ndnl:end, :) = dFdws(:, di, :, :);
            end
        else % HYSTERETIC FORCING
            Nnlf = size(m.NLTs(ni).L, 1);
            if ~isempty(m.NLTs(ni).Lf)
                Nnlf = size(m.NLTs(ni).Lf, 2);
            end
            if Nnlf>1 || size(unlt, 2)>1
                error('Not implemented for multi-DOF dependent multi-force nonlinearities yet')
            end
            % CONSTRUCT AND STORE NMAT
            switch qptype
                case {1, 2, 3, 10, 20, 30}
                    [Nmat, bpis, bpjs, Nmat_dw] = CONSTRUCTNMAT(ws, Nc, Nt, qptype);
                    bpjs = cellfun(@(c) unique(c), bpjs, 'UniformOutput', false);
                    wgt = Nmat(1,1);  % useful for case 1
                    wgt_dw = cellfun(@(m) full(m(1,1)), Nmat_dw);
        
                    t = linspace(0, 2*pi, Nt+1); t = t(1:Nt);
                    taus = cell(Nc, 1);
                    [taus{:}] = ndgrid(t);
                    Nst = QPTIMEINTERP(cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)'), h);
                    Nsti = Nst'*2/(Nt^Nc);
                    Nsti(1,:) = Nsti(1,:)/2;
                case {4, 5}
                    % Generate Regular Grid
                    t = linspace(0, 2*pi, Nt+1); t = t(1:Nt);
                    taus = cell(Nc, 1);
                    [taus{:}] = ndgrid(t);
                    taus_dw = repmat(taus, 1, Nc);
            
                    % Stagger the Grid
                    for ic=1:Nc-1
                        taus_dw{ic,ic} = taus{end}/ws(end);
                        for jc=setdiff(1:Nc-1,ic)
                            taus_dw{ic,jc} = (zeros(repmat(Nt, 1, Nc)));
                        end
                        taus_dw{ic, Nc} = -ws(ic)/ws(end)^2*taus{end};
        
                        taus{ic} = taus{ic} + ws(ic)/ws(end)*taus{end};
                    end
                    for jc=1:Nc
                        taus_dw{Nc, jc} = (zeros(repmat(Nt, 1, Nc)));
                    end
                    tstpts = cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)');  % All the time coordinates
                    tstpts_dw = cell2mat(permute(cellfun(@(c) c(:), taus_dw, 'UniformOutput', false)', [3, 2, 1]));
        
                    [Nst, Nst_dw] = QPTIMEINTERP(tstpts, h, tstpts_dw);
                    Nsti = Nst'*2/(Nt^Nc);
                    Nsti(1,:) = Nsti(1,:)/2;
        
                    Nsti_dw = permute(Nst_dw, [2 1 3])*2/(Nt^Nc);
                    Nsti_dw(1,:,:) = Nsti_dw(1,:,:)/2;
            
                    [Nmat, Nmat_dw] = CONSTRUCTNMAT(ws, Nc, Nt, 4);
                    [Nmati, Nmati_dw] = CONSTRUCTNMAT(ws, Nc, Nt, 5);
            end            

            % EVALUATE NONLINEARITY
            cst = Nst;
            ut = Nst*Uh;            
            % Trivial initial guess - works perfectly for all cases.
            ft = zeros(Nt^Nc, 1);
            dfdai = zeros(Nt^Nc, Nhc);
            dfdw = zeros(Nt^Nc, Nc);
            switch qptype
              case 1  % Approach 1 (implicit interpolation).
                x = ft(bpjs{1});  % First layer dependents
                dxdai = dfdai(bpjs{1}, :);
                dxdw = dfdw(bpjs{1}, :);

                its = 0;
                res = 0;
                while its==0 || rms(res)/muN>tol 
                    % Store first layer
                    ft(bpjs{1}) = x; 
                    dfdai(bpjs{1}, :) = dxdai; 
                    dfdw(bpjs{1}, :) = dxdw; 

                    dfpred = eye(length(x));  % Gradient wrt first layer dependents
                    dfprev = dfpred;
                    % March for one cycle
                    for ti=1:Nt
                        uc = ut(bpis{ti});
                        up = ut(bpjs{ti});
                        fc = ft(bpis{ti});
                        fp = ft(bpjs{ti});
                        cstc = cst(bpis{ti},:);
                        cstp = cst(bpjs{ti},:);
                        dfcdai = dfdai(bpis{ti},:);
                        dfpdai = dfdai(bpjs{ti},:);
                        dfcdw = dfdw(bpis{ti},:);
                        dfpdw = dfdw(bpjs{ti},:);
                        [fc, dfcdai, dfcdw, dfcdfp] = QPMARCHIT(fc, fp, uc, up, dfcdai, dfpdai, dfcdw, dfpdw, cstc, cstp, prm, ti);

                        ft(bpis{ti}) = fc;
                        dfdai(bpis{ti}, :) = dfcdai;
                        dfdw(bpis{ti}, :) = dfcdw;

                        tmp = dfcdfp*dfprev;
                        [~, si1, si2] = intersect(bpjs{mod(ti+1-1,Nt)+1}, bpis{ti});
                        [~, qi1, qi2] = intersect(bpjs{mod(ti+1-1,Nt)+1}, bpjs{ti});
                        
                        dfprev([si1;qi1],:) = [tmp(si2,:);dfprev(qi2,:)];
                        
                        [~, si1, si2] = intersect(bpjs{1}, bpjs{mod(ti+1-1,Nt)+1});
                        dfpred(si1,:) = dfprev(si2,:);
                    end
                    % Construct residue & jacobian
                    res = x-ft(bpjs{1});
                    dresdx = eye(length(x))-dfpred;
                    dresdai = dxdai - dfdai(bpjs{1},:);
                    dresdw = dxdw - dfdw(bpjs{1},:);

                    % Newton March update
                    x = x - dresdx\res;
                    dxdai = dxdai - dresdx\dresdai;
                    dxdw = dxdw - dresdx\dresdw;
                    % Regular march update if dresdx is noninvertible or so
                    % (can happen if fully stuck)
                    if ~isempty(find(isnan(x), 1))
                        x = ft(bpjs{1});
                        dxdai = dfdai(bpjs{1}, :);
                        dxdw = dfdw(bpjs{1}, :);
                    end
                    its = its+1;
                    if its>20
                        break;
                    end
                end
                %             fprintf('%d\n', its);
              case 2  % Approaches 2 (finite difference)
                x = ft(bpjs{1});
                dxdai = dfdai(bpjs{1}, :);
                dxdw = dfdw(bpjs{1}, :);
                its = 0;
                res = 0;
                while its==0 || rms(res)/muN>tol 
                    % Store first layer dependency
                    ft(bpjs{1}) = x;
                    dfdai(bpjs{1}, :) = dxdai;
                    dfdw(bpjs{1}, :) = dxdw;
                    
                    %                 dfprvs = zeros(Nt^Nc, length(x));
                    %                 dfprvs(bpjs{1}, :) = eye(length(x));
                    dfpred = eye(length(x));
                    dfprev = eye(length(x));
                    % March for one cycle
                    for ti=1:Nt
                        uc = ut(bpis{ti});  % "current" u
                        up = ut(bpjs{ti});  % "previous" u
                        fp = ft(bpjs{ti});  % "previous" f

                        fsp = kt*(uc-Nmat(bpis{ti},bpjs{ti})*up)+Nmat(bpis{ti},bpjs{ti})*fp;  % "stick prediction"
                        dfsp = kt*(cst(bpis{ti},:)-Nmat(bpis{ti},bpjs{ti})*cst(bpjs{ti},:)) + ...
                               Nmat(bpis{ti},bpjs{ti})*dfdai(bpjs{ti},:);
                        dfspdw = cell2mat(arrayfun(@(a) kt*(-Nmat_dw{a}(bpis{ti},bpjs{ti})*up)+Nmat_dw{a}(bpis{ti},bpjs{ti})*fp, 1:Nc, 'UniformOutput', false)) + ...
                                 Nmat(bpis{ti},bpjs{ti})*dfdw(bpjs{ti},:);
                        
                        dfsp(abs(fsp)>=muN, :) = 0;
                        dfspdw(abs(fsp)>=muN, :) = 0;
                        fsp(abs(fsp)>=muN) = muN*sign(fsp(abs(fsp)>=muN));

                        ft(bpis{ti}) = fsp;
                        dfdai(bpis{ti}, :) = dfsp;
                        dfdw(bpis{ti}, :) = dfspdw;

                        %                     tmp = Nmat(bpis{ti}, bpjs{ti});
                        %                     tmp(abs(fsp)>=muN, :) = 0;
                        %                     dfprvs(bpis{ti}, :) = tmp*dfprvs(bpjs{ti}, :);

                        tmp = Nmat(bpis{ti}, bpjs{ti});
                        tmp(abs(fsp)>=muN, :) = 0;
                        tmp = tmp*dfprev;

                        [~, si1, si2] = intersect(bpjs{1}, bpis{ti});
                        dfpred(si1, :) = tmp(si2, :);

                        [~, si1, si2] = intersect(bpjs{mod(ti+1-1,Nt)+1}, bpis{ti});
                        dfprev = zeros(length(bpjs{mod(ti+1-1,Nt)+1}), length(bpjs{1}));
                        dfprev(si1, :) = tmp(si2, :);
                    end
                    % Construct residue and jacobian
                    res = x - ft(bpjs{1});
                    dresdx = eye(length(x)) - dfpred;
                    %                 dresdx = eye(length(x)) - dfprvs(bpjs{1},:);
                    dresdai = dxdai - dfdai(bpjs{1}, :);
                    dresdw = dxdw - dfdw(bpjs{1}, :);

                    % Newton March update
                    x = x - dresdx\res;
                    dxdai = dxdai - dresdx\dresdai;
                    dxdw = dxdw - dresdx\dresdw;
                    % Regular march update if dresdx is noninvertible or so
                    % (can happen if fully stuck)
                    if ~isempty(find(isnan(x), 1))
                        x = ft(bpjs{1});
                        dxdai = dfdai(bpjs{1}, :);
                        dxdw = dfdw(bpjs{1}, :);
                    end
                    its = its+1;
                    if its>20
                        break;
                    end
                end
                %             fprintf('%d\n', its);
              case 3  % Approaches 3 (front marching)
                x = ft(bpjs{1});
                dxdai = dfdai(bpjs{1}, :);
                dxdw = dfdw(bpjs{1}, :);
                its = 0;
                res = 0;
                while its==0 || rms(res)/muN>tol 
                    % Store first layer dependency
                    ft(bpjs{1}) = x;
                    dfdai(bpjs{1}, :) = dxdai;
                    dfdw(bpjs{1}, :) = dxdw;
                    
                    %                 dfprvs = zeros(Nt^Nc, length(x));
                    %                 dfprvs(bpjs{1}, :) = eye(length(x));
                    dfpred = eye(length(x));
                    dfprev = eye(length(x));
                    dfprevi = zeros(length(bpis{Nt}), length(x));
                    [~, si1, si2] = intersect(bpis{Nt}, bpjs{1});
                    dfprevi(si1, :) = dfprev(si2, :);
                    % March for one cycle
                    for ti=1:Nt
                        uc = ut(bpis{ti});  % "current" u
                        up = ut(bpjs{ti});  % "previous" u
                        fp = ft(bpjs{ti});  % "previous" f

                        fsp = kt*(uc-Nmat(bpis{ti},bpjs{ti})*up)+Nmat(bpis{ti},bpjs{ti})*fp;  % "stick prediction"
                        dfsp = kt*(cst(bpis{ti},:)-Nmat(bpis{ti},bpjs{ti})*cst(bpjs{ti},:)) + ...
                               Nmat(bpis{ti},bpjs{ti})*dfdai(bpjs{ti},:);
                        dfspdw = cell2mat(arrayfun(@(a) kt*(-Nmat_dw{a}(bpis{ti},bpjs{ti})*up)+Nmat_dw{a}(bpis{ti},bpjs{ti})*fp, 1:Nc, 'UniformOutput', false)) + ...
                                 Nmat(bpis{ti},bpjs{ti})*dfdw(bpjs{ti},:);
                        
                        dfsp(abs(fsp)>=muN, :) = 0;
                        dfspdw(abs(fsp)>=muN, :) = 0;
                        fsp(abs(fsp)>=muN) = muN*sign(fsp(abs(fsp)>=muN));

                        ft(bpis{ti}) = fsp;
                        dfdai(bpis{ti}, :) = dfsp;
                        dfdw(bpis{ti}, :) = dfspdw;

                        %                     tmp = Nmat(bpis{ti}, bpjs{ti});
                        %                     tmp(abs(fsp)>=muN, :) = 0;
                        %                     dfprvs(bpis{ti}, :) = tmp*dfprvs(bpjs{ti}, :);

                        tmp = Nmat(bpis{ti}, bpjs{ti});
                        tmp(abs(fsp)>=muN, :) = 0;
                        tmp = tmp*dfprev; % der(i, j1)

                        [~, si1, si2] = intersect(bpjs{1}, bpis{ti});
                        dfpred(si1, :) = tmp(si2, :);

                        [~, si1, si2] = intersect(bpjs{1}, bpjs{mod(ti+1-1,Nt)+1});
                        [~, qi1, qi2] = intersect(bpis{ti}, bpjs{mod(ti+1-1,Nt)+1});
                        [~, ri1, ri2] = intersect(bpis{mod(ti-1-1,Nt)+1}, bpjs{mod(ti+1-1,Nt)+1});

                        dfprev = zeros(length(bpjs{mod(ti+1-1,Nt)+1}), length(bpjs{1}));
                        dfprev(ri2, :) = dfprevi(ri1, :);
                        dfprev(si2, :) = dfpred(si1, :);
                        dfprev(qi2, :) = tmp(qi1, :);

                        dfprevi = tmp;
                    end
                    % Construct residue and jacobian
                    res = x - ft(bpjs{1});
                    dresdx = eye(length(x)) - dfpred;
                    %                 dresdx = eye(length(x)) - dfprvs(bpjs{1},:);
                    dresdai = dxdai - dfdai(bpjs{1}, :);
                    dresdw = dxdw - dfdw(bpjs{1}, :);

                    % Newton March update
                    x = x - dresdx\res;
                    dxdai = dxdai - dresdx\dresdai;
                    dxdw = dxdw - dresdx\dresdw;
                    % Regular march update if dresdx is noninvertible or so
                    % (can happen if fully stuck)
                    if ~isempty(find(isnan(x), 1))
                        x = ft(bpjs{1});
                        dxdai = dfdai(bpjs{1}, :);
                        dxdw = dfdw(bpjs{1}, :);
                    end
                    its = its+1;
                    if its>20
                        break;
                    end
                end
                %             fprintf('%d\n', its);
              case {4, 5}  % March on Stagger-Stretched Grid
                ut_dw = cell2mat(arrayfun(@(a) Nst_dw(:,:,a)*Uh,1:Nc,'UniformOutput',false));
                
                ipr = (Nt-1)*(Nt^(Nc-1))+(1:Nt^(Nc-1));
                x = ft(ipr);
                dxdai = dfdai(ipr, :);
                dxdw = dfdw(ipr, :);
                
                its = 0;
                res = 0;
                while its==0 || rms(res)/muN>tol
                    % Store first layer dependents
                    ipr = (Nt-1)*(Nt^(Nc-1))+(1:Nt^(Nc-1));
                    icr = 1:Nt^(Nc-1);
                    ft(ipr) = x;  % Original first layer
                    dfdai(ipr, :) = dxdai;
                    dfdw(ipr, :) = dxdw;

                    % Update "initial plane"                    
                    fsp = kt*(Nmat*ut(icr) - ut(ipr)) + ft(ipr);
                    dfsp = kt*(Nmat*Nst(icr,:)-Nst(ipr,:)) + dfdai(ipr,:);
                    dfspdw = cell2mat(arrayfun(@(a) kt*(Nmat_dw{a}*ut(icr)), 1:Nc, 'UniformOutput', false)) + ...
                             kt*(Nmat*ut_dw(icr,:)-ut_dw(ipr,:)) + dfdw(ipr,:);
                    
                    dfsp(abs(fsp)>=muN, :) = 0;
                    dfspdw(abs(fsp)>=muN, :) = 0;
                    fsp(abs(fsp)>=muN) = muN*sign(fsp(abs(fsp)>=muN));
                    
                    ft(icr) = Nmati*fsp;
                    dfdai(icr,:) = Nmati*dfsp;
                    dfdw(icr,:) = Nmati*dfspdw + ...
                        cell2mat(arrayfun(@(a) Nmati_dw{a}*fsp, 1:Nc, 'UniformOutput', false));
                    
                    dfpred = Nmati*diag(abs(fsp)<muN);
                    % March over rest of the points for one cycle
                    for ti=2:Nt
                        ipr = icr;  % prev. pt indices
                        icr = (ti-1)*(Nt^(Nc-1))+(1:Nt^(Nc-1));  % current pt indices

                        fsp = kt*(ut(icr) - ut(ipr)) + ft(ipr);
                        dfsp = kt*(Nst(icr,:)-Nst(ipr,:)) + dfdai(ipr,:);
                        dfspdw = kt*(ut_dw(icr,:) - ut_dw(ipr,:)) + dfdw(ipr,:);

                        dfsp(abs(fsp)>=muN,:) = 0;
                        dfspdw(abs(fsp)>=muN,:) = 0;
                        fsp(abs(fsp)>=muN) = muN*sign(fsp(abs(fsp)>=muN));

                        ft(icr) = fsp;
                        dfdai(icr,:) = dfsp;
                        dfdw(icr,:) = dfspdw;

                        dfpred = diag(abs(fsp)<muN)*dfpred;
                    end
                    % Construct residue
                    ipr = icr;
                    icr = 1:Nt^(Nc-1);

                    res = x - ft(ipr);
                    dresdx = eye(length(x)) - dfpred;     
                    dresdai = dxdai - dfdai(ipr, :);
                    dresdw = dxdw - dfdw(ipr, :);

                    % Newton update march
                    x = x - dresdx\res;
                    dxdai = dxdai - dresdx\dresdai;
                    dxdw = dxdw - dresdx\dresdw;
                    % Regular march update if dresdx is noninvertible or so
                    % (can happen if fully stuck)
                    if ~isempty(find(isnan(x),1))
                        x = ft(ipr);
                        dxdai = dfdai(ipr, :);
                        dxdw = dfdw(ipr, :);
                    end
                    its = its+1;
                    if its>20
                        break;
                    end
                end
              case {10,20,30}  % Solve using fsolve. Versions of case 1, 2.
                cst = QPAFT([eye(Nhc) Uh], h, Nt, 'f2t');
                ut = cst(:, Nhc+1:end);
                cst = cst(:, 1:Nhc);

                opt = optimoptions('fsolve', 'specifyobjectivegradient', true, 'Display', 'iter');
                ft = ones(Nt^Nc, 1);  % Better Initial guess?
                ft = fsolve(@(ft) QPMARCHRESFUN(ft, ut, prm), ft, opt);

                % Get Jacobians
                [~, dresdf, dresdu, dresdw] = QPMARCHRESFUN(ft, ut, prm);
                dfdai = (-dresdf\dresdu)*cst;
                dfdw = -dresdf\dresdw;
            end
            
            
            switch m.NLTs(ni).qptype
                case {1,2}  % Option 1: Solve using NSOLVE or fsolve
                    % Construct Nmat
                    ft = ones(Nt^Nc, Nnlf);
                    
                    opt = struct('Display', false);
        %             tic
                    ft = NSOLVE(@(ft) QPMARCHRESFUN(ft, unlt, ...
                        @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
                        Nmat), ft, opt);
        %             toc
        %             opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'off');
        %             ft = fsolve(@(ft) QPMARCHRESFUN(ft, unlt, ...
        %                 @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
        %                 Nmat), ft, opt);
        
                    % Get Jacobians
                    [~, dresdf, dresdu] = QPMARCHRESFUN(ft, unlt, ...
                        @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
                        Nmat);
                    dfdut = -dresdf\dresdu;
                    
                    F = QPFOURIERCOEFF(ft, h);
                    J = QPFOURIERCOEFF(dfdut*cst, h);
                case 3  % Option 2: Sequential Solution: Only possible for Nmtype=3
                    [Nmat, bpis, bpjs] = CONSTRUCTNMAT(Ws, Nc, Nt, m.NLTs(ni).qptype);
                    bpjs = cellfun(@(c) unique(c), bpjs, 'UniformOutput', false);
        
                    its = 0;
                    ft = zeros(Nt^Nc, Nnlf);
                    dfdai = zeros(Nt^Nc, Nhc);
                    fprev = ft(bpis{1},:);
        %             tic
                    while its==0 || max(abs(ft(bpis{1})-fprev))>tol 
                        fprev = ft(bpis{1},:);
                        for ti=1:Nt
                            % Only force estimation
                            ft(bpis{ti}) = m.NLTs(ni).func(0, unlt(bpis{ti}), 0, 0, Nmat(bpis{ti},bpjs{ti})*unlt(bpjs{ti}), Nmat(bpis{ti},bpjs{ti})*ft(bpjs{ti}));
        
        %                     % Force & Jacobian estimation - Naive version
        %                     [ft(bpis{ti}), ~, dfdfp, dfdu, dfdup] = m.NLTs(ni).func(0, unlt(bpis{ti}), 0, 0, Nmat(bpis{ti},:)*unlt, Nmat(bpis{ti},:)*ft);
        %                     dfdai(bpis{ti},:) = diag(dfdfp)*Nmat(bpis{ti},:)*dfdai + dfdu.*cst(bpis{ti},:) + diag(dfdup)*Nmat(bpis{ti},:)*cst;
        
        %                     % Force & Jacobian estimation - respectful of sparsity in Nmat
        %                     [ft(bpis{ti}), ~, dfdfp, dfdu, dfdup] = m.NLTs(ni).func(0, unlt(bpis{ti}), 0, 0, Nmat(bpis{ti},bpjs{ti})*unlt(bpjs{ti}), Nmat(bpis{ti},bpjs{ti})*ft(bpjs{ti}));
        %                     dfdai(bpis{ti},:) = diag(dfdfp)*Nmat(bpis{ti},bpjs{ti})*dfdai(bpjs{ti},:) + dfdu.*cst(bpis{ti},:) + diag(dfdup)*Nmat(bpis{ti},bpjs{ti})*cst(bpjs{ti},:);
                        end
                        its = its+1;
                    end
                    % Iterate and just estimate the Jacobians
                    for ti=1:Nt  % Run the iterations once more
                        [ft(bpis{ti}), ~, dfdfp, dfdu, dfdup] = m.NLTs(ni).func(0, unlt(bpis{ti}), 0, 0, Nmat(bpis{ti},bpjs{ti})*unlt(bpjs{ti}), Nmat(bpis{ti},bpjs{ti})*ft(bpjs{ti}));
                        dfdai(bpis{ti},:) = diag(dfdfp)*Nmat(bpis{ti},bpjs{ti})*dfdai(bpjs{ti},:) + dfdu.*cst(bpis{ti},:) + diag(dfdup)*Nmat(bpis{ti},bpjs{ti})*cst(bpjs{ti},:);
                    end
%                     fprintf('%d\n',its)
                    F = QPFOURIERCOEFF(ft, h);
                    J = QPFOURIERCOEFF(dfdai, h);
        %             toc
            end
            
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
