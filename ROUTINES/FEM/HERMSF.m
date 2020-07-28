function [Ns, dNs] = HERMSF(xi, Le)
    Ns = [2*(xi-1).^2.*(xi+2), Le*(xi-1).^2.*(xi+1),...
        -2*(xi-2).*(xi+1).^2, Le*(xi-1).*(xi+1).^2]/8;
    
    dNs = [(3*(xi-1)*(xi+1))/4, (Le*(xi-1)*(3*xi+1))/8, ...
        -(3*(xi-1)*(xi+1))/4, (Le*(xi+1)*(3*xi-1))/8];
    
    % at 0: [0.5, Le/8, 0.5, -Le/8]
end