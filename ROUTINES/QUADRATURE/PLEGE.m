function [p, intp] = PLEGE(n, x)
    if length(n)==1
        if n<0
            p = zeros(size(x));
        elseif n==0
            p = ones(size(x));
        else
            p = ((2*n-1)*x.*PLEGE(n-1, x) - (n-1).*PLEGE(n-2, x))/n;
        end
        
        if n<0
            intp = 0;
        else
            intp = 2/(2*n+1);
        end
    else
        h = zeros(length(x), length(n));
        inth = zeros(length(n), 1);
        if strcmp(class(x), 'sym')
            h = sym(h);
        end
        for in=1:length(n)
            [h(:, in), inth(in)] = PLEGE(n(in), x);
        end
    end
end