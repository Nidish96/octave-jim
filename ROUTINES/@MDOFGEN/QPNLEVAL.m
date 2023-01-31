function [FNL, dFNL, dFNLdw] = QPNLEVAL(m, Uws, h, Nt, tol)
%QPNLEVAL evaluates nonlinear force in QP-frequency domain
%
%   USAGE:
%       [FNL, dFNL, dFNLdw] = QPNLEVAL(m, Uws, h, Nt, tol);
%   INPUTS:
%       Uws         : (Nd*Nhc+Nc,1)
%       h           : (Nh, Nc)
%       Nt          : (1,1)
%       tol         : (1,1)
%   OUTPUTS:
%       FNL         : (Nd*Nhc,1)
%       dFNL        : (Nd*Nhc, Nd*Nhc)
%       dFNLdw      : (Nd*Nhc, Nc)


    Nc = size(h,2);  % Number of components
    Nh = size(h,1);  % Number of harmonics
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));  % Number of harmonic coefficients
    
    ws = Uws(end-Nc+1:end);  % List of frequencies

    tau = linspace(0, 2*pi, Nt+1); tau(end) = [];
    taus = cell(Nc, 1);
    [taus{:}] = ndgrid(tau);
    t = zeros(repmat(Nt, 1, Nc));
    for ic=1:Nc
        t = t+taus{ic}*ws(ic);
    end    
    [D1, dD1dws] = QPHARMONICSTIFFNESS(0, 1, 0, ws, h);  % time derivative matrix
    
    cst = QPTIMETRANS(eye(Nhc), h, Nt);  % basis functions
    sct = QPTIMETRANS(D1, h, Nt);  % basis function derivatives
    dsct_ws = cell2mat(permute(arrayfun(@(a) QPTIMETRANS(dD1dws(:,:,a), h, Nt), 1:Nc, 'UniformOutput', false), [1 3 2]));

    FNL = zeros(m.Ndofs*Nhc, 1);
    dFNL = zeros(m.Ndofs*Nhc);
    dFNLdw = zeros(m.Ndofs*Nhc, Nc);
    for ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(Uws(1:end-Nc), m.Ndofs, Nhc))';  % Nhc x Ndnl
        Ndnl = size(m.NLTs(ni).L, 1);
        
        if mod(m.NLTs(ni).type, 2)==0  % INSTANTANEOUS FORCING
            unlt   = QPTIMETRANS(Unl, h, Nt);
            unldot = QPTIMETRANS(D1*Unl, h, Nt);

            [ft, dfdu, dfdud] = m.NLTs(ni).func(taus, unlt, unldot);
        
            F = QPFOURIERCOEFF(ft, h);
            J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
            dFdU = reshape(QPFOURIERCOEFF(reshape(dfdu.*permute(cst, [1, 3, 2]), Nt^Nc, Ndnl*Nhc), h) + ...
                QPFOURIERCOEFF(reshape(dfdud.*permute(sct, [1, 3, 2]), Nt^Nc, Ndnl*Nhc), h), ...
                Nhc, Ndnl, Nhc);
            % Hasn't been generalized to mutliple DOFs yet
            dFdws = reshape(cell2mat(arrayfun(@(a) QPFOURIERCOEFF(reshape(dfdud.*permute(dsct_ws(:,:,a), [1, 3, 2]), Nt^Nc, Ndnl*Nhc),h)*Unl, 1:Nc, 'UniformOutput', false)), ...
                Nhc, Ndnl, Nc);
            
            Jw = zeros(length(F), Nc);
            for di=1:Ndnl
                J(di:Ndnl:end, di:Ndnl:end) = dFdU(:, di, :);
                Jw(di:Ndnl:end, :) = dFdws(:, di, :);
            end
        else % HYSTERETIC FORCING
            Nnlf = size(m.NLTs(ni).L, 1);
            if ~isempty(m.NLTs(ni).Lf)
                Nnlf = size(m.NLTs(ni).Lf, 2);
            end
            if Nnlf>1 || size(Unl, 2)>1
                error('Not implemented for multi-DOF dependent multi-force nonlinearities yet')
            end

            %% Preprocessing: Construct and Store Nmat
            switch m.NLTs(ni).qptype
              case {1, 2, 3, 10, 20, 30}
                [Nmat, bpis, bpjs, Nmat_dw] = CONSTRUCTNMAT(ws(:)', Nc, Nt, m.NLTs(ni).qptype);
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
                
                [Nmat, Nmat_dw] = CONSTRUCTNMAT(ws(:)', Nc, Nt, 4);
                [Nmati, Nmati_dw] = CONSTRUCTNMAT(ws(:)', Nc, Nt, 5);
            end

            %% Evaluate Nonlinearity
            cst = Nst;
            ut = Nst*Unl;

            % Trivial initial guess - works perfectly for all cases.
            ft = zeros(Nt^Nc, 1);
            dfdai = zeros(Nt^Nc, Nhc);
            dfdw = zeros(Nt^Nc, Nc);            
            
            switch m.NLTs(ni).qptype
              case {4, 5}  % March on Stagger-Stretched Grid
                ut_dw = cell2mat(arrayfun(@(a) Nst_dw(:,:,a)*Unl,1:Nc,'UniformOutput',false));
                
                ipr = (Nt-1)*(Nt^(Nc-1))+(1:Nt^(Nc-1));
                x = ft(ipr);
                dxdai = dfdai(ipr, :);
                dxdw = dfdw(ipr, :);
                
                its = 0;
                res = 0;
                while its==0 || rms(res)>tol
                    % Store first layer dependents
                    ipr = (Nt-1)*(Nt^(Nc-1))+(1:Nt^(Nc-1));
                    icr = 1:Nt^(Nc-1);
                    ft(ipr) = x;  % Original first layer
                    dfdai(ipr, :) = dxdai;
                    dfdw(ipr, :) = dxdw;

                    % Update "initial plane"
                    [fc, ~, dfdfp, dfduc, dfdup] = m.NLTs(ni).func(1, Nmat*ut(icr), 0, 0, ...
                                                                   ut(ipr), ft(ipr));
                    ft(icr) = Nmati*fc;
                    dfdai(icr, :) = Nmati*(dfduc.*(Nmat*Nst(icr,:)) + ...
                                           dfdup.*Nst(ipr,:) + ...
                                           dfdfp.*dfdai(ipr, :));
                    dfdw(icr, :) = cell2mat(arrayfun(@(a) Nmati*(dfduc.*Nmat_dw{a}*ut(icr)) + ...
                                                     Nmati_dw{a}*fc, 1:Nc, ...
                                                     'Uniformoutput', false)) + ...
                                                     Nmati*(dfdfp.*dfdw(ipr,:) + ...
                                                     dfduc.*Nmat*ut_dw(icr,:) + ...
                                                     dfdup.*ut_dw(ipr,:));
                    if any(dfdfp==0)
                        mrflag=1;
                    else
                        mrflag=0;
                    end
                    dfpred = Nmati*diag(dfdfp);
                    % March over rest of the points for one cycle
                    for ti=2:Nt
                        ipr = icr;  % prev. pt indices
                        icr = (ti-1)*(Nt^(Nc-1))+(1:Nt^(Nc-1));  % current pt indices

                        [fc, ~, dfdfp, dfduc, dfdup] = m.NLTs(ni).func(ti, ut(icr), 0, 0, ...
                                                                       ut(ipr), ft(ipr));
                        ft(icr) = fc;
                        dfdai(icr, :) = dfduc.*Nst(icr,:) + dfdup.*Nst(ipr,:) + dfdfp.*dfdai(ipr,:);
                        dfdw(icr, :) = dfduc.*ut_dw(icr, :) + dfdup.*ut_dw(ipr, :) + ...
                            dfdfp.*dfdw(ipr, :);

                        if any(dfdfp==0)  % slip-like memory resetting event
                            mrflag = 1;
                        end

                        dfpred = dfdfp.*dfpred;
                    end
                    if mrflag==0  % No slip-like events. Fully stuck.
                        break;
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
                        warning('no slip');
                        x = ft(ipr);
                        dxdai = dfdai(ipr, :);
                        dxdw = dfdw(ipr, :);
                    end
                    its = its+1;
                    if its>20
                        break;
                    end
                end
              otherwise
                error('Other cases not implemented here yet');
            end
            
            F = Nsti*ft;
            J = Nsti*dfdai;
            Jw = Nsti*dfdw + ...
                 cell2mat(arrayfun(@(a) Nsti_dw(:, :, a)*ft, 1:Nc, 'Uniformoutput', false));
        end
        
        if m.NLTs(ni).type<=5  % Self adjoint forcing
            FNL  = FNL + reshape(m.NLTs(ni).L'*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).L')*J*kron(eye(Nhc), m.NLTs(ni).L);
            dFNLdw = dFNLdw + kron(eye(Nhc), m.NLTs(ni).L')*Jw;
        else  % Non-self adjoint forcing
            FNL  = FNL + reshape(m.NLTs(ni).Lf*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).Lf)*J*kron(eye(Nhc), m.NLTs(ni).L);
            dFNLdw = dFNLdw + kron(eye(Nhc), m.NLTs(ni).Lf)*Jw;
        end
    end
end
