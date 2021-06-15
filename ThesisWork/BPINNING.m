function [f, dfdu, dfdud] = BPINNING(t, u, ud, ks, clrs)
%BPINNING returns the pinned bolt forces
%
%   USAGE :
%       [f, dfdu, dfdud] = BPINNING(t, u, ud, ks, clrs);
%   INPUTS :
%       t   :
%       u   :
%       ud  :
%       ks  :
%       clrs:
    Np = length(u)/2;
    
    ks = ks(:);  clrs = clrs(:);
    if length(ks)==1 || length(clrs)==1
        ks   = ks*ones(Np,1);
        clrs = clrs*ones(Np,1);
    end
    
    un = sqrt(u(1:2:end).^2+u(2:2:end).^2);
    dundxy = zeros(Np, Np*2);
    for i=1:Np
        dundxy(i, (i-1)*2+(1:2)) = u((i-1)*2+(1:2))/un(i);
    end
    
    fn = ks.*max(un-clrs, 0);
    dfndu = diag(ks.*(un>clrs))*dundxy;
    
    f = kron(fn./un, [1;1]).*u;
    dfdu = diag(kron(fn./un, [1;1])) + ...
        (kron(dfndu./un, [1;1]) - ...
        kron(fn./un.^2.*dundxy, [1;1])).*u;
    dfdud = dfdu*0;
    
%     if any(f~=0)
%         keyboard
%     end
end