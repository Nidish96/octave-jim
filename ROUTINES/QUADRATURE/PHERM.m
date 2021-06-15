function [h, inth] = PHERM(n, x)
    if length(n)==1
        if n<0
            h = zeros(size(x));
        elseif n==0
            h = ones(size(x));
        else
            h = x.*PHERM(n-1, x) - (n-1)*PHERM(n-2, x);
        end
        
        if n<0
            inth = 0;
        else
            inth = factorial(n);
        end
    else
        h = zeros(length(x), length(n));
        inth = zeros(length(n), 1);
        if strcmp(class(x), 'sym')
            h = sym(h);
        end
        for in=1:length(n)
            [h(:, in), inth(in)] = PHERM(n(in), x);
        end        
    end
end