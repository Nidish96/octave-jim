function [h, inth] = PLAGU(n, x)
    if length(n)==1
        if n<0
            h = zeros(size(x));
        elseif n==0
            h = ones(size(x));
        else
            h = ((2*n-1-x).*PLAGU(n-1, x) - (n-1)*PLAGU(n-2, x))/n;
        end

        if n<0
            inth = 0;
        else
            inth = 1;
        end
    else % Multiple polynomials required            
        h = zeros(length(x), length(n));
        inth = zeros(length(n), 1);
        if strcmp(class(x), 'sym')
            h = sym(h);
        end
        for in=1:length(n)
            [h(:, in), inth(in)] = PLAGU(n(in), x);
        end
    end
end