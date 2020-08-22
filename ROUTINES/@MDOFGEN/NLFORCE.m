function [F, dFdU, dFdUd, m] = NLFORCE(m, t, U, Ud, tp, varargin)
    
    if length(varargin)==1
        init = varargin{1};  % 1 if this is the initalization run, 0 if not
    else
        init = 0;
    end
        

    F = zeros(m.Ndofs, 1);
    dFdU = zeros(m.Ndofs);
    dFdUd = zeros(m.Ndofs);
    
    for ni=1:length(m.NLTs)
        if mod(m.NLTs(ni).type,2)==0  % Inst. force
            [f, dfdu, dfdud] = m.NLTs(ni).func(t, m.NLTs(ni).L*U, m.NLTs(ni).L*Ud);
        else
            if init
                [f, dfdu] = m.NLTs(ni).func(t, m.NLTs(ni).L*U);
            else
                [f, dfdu] = m.NLTs(ni).func(t, m.NLTs(ni).L*U, 0, tp, m.NLTs(ni).up, ...
                    m.NLTs(ni).fp, m.NLTs(ni).dfdup);
            end
            dfdud = zeros(size(dfdu));

            m.NLTs(ni).up = m.NLTs(ni).L*U;
            m.NLTs(ni).fp = f;
            m.NLTs(ni).dfdup = dfdu;
        end

        if m.NLTs(ni).type<=5
            F = F + m.NLTs(ni).L'*f;
            dFdU = dFdU + m.NLTs(ni).L'*dfdu*m.NLTs(ni).L;
            dFdUd = dFdUd + m.NLTs(ni).L'*dfdud*m.NLTs(ni).L;
        else
            F = F + m.NLTs(ni).Lf*f;
            dFdU = dFdU + m.NLTs(ni).Lf*dfdu*m.NLTs(ni).L;
            dFdUd = dFdUd + m.NLTs(ni).Lf*dfdud*m.NLTs(ni).L;
        end
    end
end