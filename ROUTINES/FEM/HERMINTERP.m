function [Uq, Nint, dNint, qis, eis, xis] = HERMINTERP(X, U, Xq, varargin)
%HERMINTERP interpolates the set of values and gradients using Hermite
%Polynomials
%  
%  USAGE:
%       Xq = HERMINTERP(X, U, Xq);
%  INPUTS:
%       X   : (N,1)
%       U   : (2N,1)
%       Xq  : (Nq,1)
%       npar: number of threads to run for the interpolation function
%  OUTPUTS:
%       Uq  : (Nq,1)

    if nargin==4
        npar = varargin{1};
    else
        npar = 6;
    end
    
    U = U(:);
    
    Uq = zeros(length(Xq), 1);
    
    X = X(:);
    Xq = Xq(:);
    
    Les = diff(X);
    Ne  = length(Les);
    
    [qis,eis] = find((X(1:end-1)'-Xq).*(X(2:end)'-Xq)<=0);
    [qis, si] = sort(qis);    eis = eis(si);
    [qis, si] = unique(qis);  eis = eis(si);
    
    if sum(qis(:)'-(1:length(qis)))~=0
        keyboard
        error('no way')
    end
    
    xis = zeros(size(qis));
    Nint = zeros(length(Xq), length(X)*2);
    dNint = zeros(length(Xq), length(X)*2);
    
    opt = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
    % Serial
%     for qi=1:length(qis)
%         e = eis(qi);
%         xis(qi) = fsolve(@(xi) RESFUN(xi, [X(e);1.0;X(e+1);1.0], Xq(qi), Les(e)), 0, opt);
%         Uq(qi) = HERMSF(xis(qi), Les(e))*U((e-1)*2+1:(e+1)*2);
%         
%         [Nint(qi, (e-1)*2+1:(e+1)*2), dNint(qi, (e-1)*2+1:(e+1)*2)] = HERMSF(xis(qi), Les(e));
%         dNint(qi, (e-1)*2+1:(e+1)*2) = dNint(qi, (e-1)*2+1:(e+1)*2)*2/Les(e);
%     end

    % Parallel
    Ues = [U(1:2:end-3) U(2:2:end-2) U(3:2:end-1) U(4:2:end)];
    Xes = [X(1:end-1) ones(Ne,1) X(2:end) ones(Ne,1)];
    Les = diff(X);
    
    Xes = Xes(eis, :);
    Ues = Ues(eis, :);
    Les = Les(eis);
    parfor (qi=1:length(qis), npar)
        e = eis(qi);
        xis(qi) = fsolve(@(xi) RESFUN(xi, Xes(qi,:)', Xq(qi), Les(qi)),...
            0, opt);
        [N, dN] = HERMSF(xis(qi), Les(qi));
        
        Uq(qi) = Ues(qi, :)*N';
        
        Nint(qi, :) = [zeros(1,(e-1)*2), N, zeros(1,(Ne+1)*2-(e+1)*2)];
        dNint(qi, :) = [zeros(1,(e-1)*2), dN*2/Les(qi), zeros(1,(Ne+1)*2-(e+1)*2)];
    end    
    Nint = sparse(Nint);
    dNint = sparse(dNint);
end

function [R, dRdxi] = RESFUN(xi, Xe, Xq, Le)
    [Ns, dNs] = HERMSF(xi, Le);
    
    R = Ns*Xe-Xq;
    dRdxi = dNs*Xe;
end