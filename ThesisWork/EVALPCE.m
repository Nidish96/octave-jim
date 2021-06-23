function [y] = EVALPCE(x, Ycofs, IJs, polfuns)
%EVALPCE evaluates the given PCE at required points x

    Psi = 1;  % Bases
    for j=1:size(IJs, 2)
        psi = polfuns{j}(IJs(:,j), x(:, j));
        Psi = Psi.*psi;
    end
    y = Ycofs*Psi';
end