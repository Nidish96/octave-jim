function [Uq, Nint, dNint, qis, eis, xis] = HERMINTERP(X, U, Xq)
%HERMINTERP interpolates the set of values and gradients using Hermite
%Polynomials
%  
%  USAGE:
%       Xq = HERMINTERP(X, U, Xq);
%  INPUTS:
%       X   : (N,1)
%       U   : (2N,1)
%       Xq  : (Nq,1)
%  OUTPUTS:
%       Uq  : (Nq,1)

    U = U(:);
    
    Uq = zeros(length(Xq), 1);
    
    X = X(:);
    Xq = Xq(:);
    
    [qis,eis] = find((X(1:end-1)'-Xq).*(X(2:end)'-Xq)<=0);
    [qis, si] = sort(qis);    eis = eis(si);
    [qis, si] = unique(qis);  eis = eis(si);
    
    xis = zeros(size(qis));
    Nint = zeros(length(Xq), length(X)*2);
    dNint = zeros(length(Xq), length(X)*2);
    
    opt = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
    for i=1:length(qis)
        qi = qis(i);
        e = eis(i);
        Le = X(e+1)-X(e);
        xis(i) = fsolve(@(xi) RESFUN(xi, [X(e);1.0;X(e+1);1.0], Xq(qi), Le), 0, opt);
        Uq(qi) = HERMSF(xis(i), Le)*U((e-1)*2+1:(e+1)*2);
        
        [Nint(qi, (e-1)*2+1:(e+1)*2), dNint(qi, (e-1)*2+1:(e+1)*2)] = HERMSF(xis(i), Le);
        dNint(qi, (e-1)*2+1:(e+1)*2) = dNint(qi, (e-1)*2+1:(e+1)*2)*2/Le;
    end
    Nint = sparse(Nint);
    dNint = sparse(dNint);
end

function [R, dRdxi] = RESFUN(xi, Xe, Xq, Le)
    [Ns, dNs] = HERMSF(xi, Le);
    
    R = Ns*Xe-Xq;
    dRdxi = dNs*Xe;
end