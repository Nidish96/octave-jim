function [FY, fY] = POLYCDF(gp, y, FX, fX, lub, varargin)
%     x = sym('x');
%     g = poly2sym(gp, x);
%     xs = double(solve(g==y, x));
%     dydxs = double(subs(diff(g, x), {x}, {xs}));

    gs = gp; gs(end) = gs(end)-y;
    rs = roots(gs); 
    xs = rs(imag(rs)==0 & rs>=lub(1) & rs<=lub(2));
    dydxs = polyval(polyder(gp), xs);
    
    % Extremal Points
    es = [0; roots(polyder(gp))];
    es = es(imag(es)==0 & es>=lub(1) & es<=lub(2));
    exm = max(es);
    d2dx2 = polyval(polyder(polyder(gp)), exm);
    
    if length(xs)==1
        FY = FX(xs)*sign(dydxs);
        fY = fX(xs)./abs(dydxs);
    else
        FY = sum(FX(xs).*sign(dydxs));
        fY = sum(fX(xs)./abs(dydxs));
    end
    FY = (sign(d2dx2)==-1)+FY;
    
    if length(varargin)==1  % Target CDF value [0, 1]
        FY = FY-varargin{1};
    end
end