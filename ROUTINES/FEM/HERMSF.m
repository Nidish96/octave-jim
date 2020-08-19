function [Ns, dNs] = HERMSF(xi, Le)
    L1 = ones(size(Le));
    
    Ns = [2*(xi-1).^2.*(xi+2).*L1, Le.*(xi-1).^2.*(xi+1),...
        -2*(xi-2).*(xi+1).^2.*L1, Le.*(xi-1).*(xi+1).^2]/8;
    
    dNs = [(6*(xi-1).*(xi+1)).*L1, (Le.*(xi-1).*(3*xi+1)), ...
        -(6*(xi-1).*(xi+1)).*L1, (Le.*(xi+1).*(3*xi-1))]/8;
    
    % at 0: [0.5, Le/8, 0.5, -Le/8]
end