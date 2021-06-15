function [FY, fY] = SYMCDF(g, y, x, FX, fX, varargin)
    xs = double(solve(g==y, x));
    dxdys = double(subs(diff(g, x), {x}, {xs}));
    
    if length(xs)==1
        FY = FX(xs);
        fY = fX(xs)./dxdys;
    else
        FY = sum(FX(xs).*sign(dxdys));
        fY = sum(fX(xs)./abs(dxdys));
    end
    
    if length(varargin)==1  % Target CDF value [0, 1]
        FY = FY-varargin{1};
    end
end