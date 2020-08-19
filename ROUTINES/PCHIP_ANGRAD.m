function [yq, pn] = PCHIP_ANGRAD(x,y,slopes,xq)
%PCHIP  Piecewise Cubic Hermite Interpolating Polynomial.
%   YQ = PCHIP(X,Y,XQ) performs shape-preserving cubic Hermite
%   interpolation using the values Y at sample points X to find
%   interpolated values YQ at the query points XQ.
%       - X must be a vector.
%       - If Y is a vector, Y(j) is the value at X(j).
%       - If Y is a matrix or n-D array, Y(:,...,:,j) is the value at X(j).
%
%   PCHIP chooses slopes at X(j) such that YQ respects monotonicity and is
%   shape-preserving. Namely, on intervals where the Y data is monotonic,
%   so is YQ; at points where the Y data has a local extremum, so does YQ.
%
%   PP = PCHIP(X,Y) returns the piecewise polynomial form PP of the
%   interpolant. You can use PP as an input to PPVAL or UNMKPP.
%
%   Comparison of SPLINE, PCHIP, and MAKIMA:
%       - All three are a form of piecewise cubic Hermite interpolation,
%         but each function computes the slopes of YQ at X(j) differently.
%       - SPLINE chooses slopes at X(j) such that the second derivative of
%         YQ is continuous. Therefore, SPLINE is smoother and more accurate
%         if the Y data represents values of a smooth function.
%       - PCHIP has no overshoots and less oscillation than SPLINE.
%       - MAKIMA has less oscillation than SPLINE but may have overshoots.
%       - PCHIP and MAKIMA are less expensive than SPLINE to set up PP.
%       - All three are equally expensive to evaluate.
%       - SPLINE and MAKIMA generalize to n-D grids. See INTERPN.
%
%   Example: Compare SPLINE, PCHIP, and MAKIMA
%
%       x = [1 2 3 4 5 5.5 7 8 9 9.5 10];
%       y = [0 0 0 0.5 0.4 1.2 1.2 0.1 0 0.3 0.6];
%       xq = 0.75:0.05:10.25;
%       yqs = spline(x,y,xq);
%       yqp = pchip(x,y,xq);
%       yqm = makima(x,y,xq);
%
%       plot(x,y,'ko','LineWidth',2,'MarkerSize',10)
%       hold on
%       plot(xq,yqp,'LineWidth',4)
%       plot(xq,yqs,xq,yqm,'LineWidth',2)
%       legend('(x,y) data','pchip','spline','makima')
%
%   See also INTERP1, MAKIMA, SPLINE, PPVAL, MKPP, UNMKPP.

    % References:
    %   F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
    %   Interpolation", SIAM J. Numerical Analysis 17, 1980, 238-246.
    %   David Kahaner, Cleve Moler and Stephen Nash, Numerical Methods
    %   and Software, Prentice Hall, 1988.

    %   Copyright 1984-2019 The MathWorks, Inc.

    % Check and adjust input data
    [x,y,sizey] = CHCKXY(x,y);
    [~,slopes,~] = CHCKXY(x,slopes);

    % Compute slopes
    h = diff(x);
    m = prod(sizey);
    delta = diff(y,1,2)./repmat(h,m,1);
    % if size(slopes,1)~=size(y,1)
    %     slopes = slopes';
    % end

    % Compute piecewise cubic Hermite interpolant for those values and slopes
    pn = pwch(x,y,slopes,h,delta);
    pn.dim = sizey;

    % Evaluate the piecewise cubic Hermite interpolant
    yq = ppval(pn,xq);
end
