function [R, dRdU, m] = STATRESFUN(m, U, Fstat)
%STATRESFUN returns the static residue & Jacobian
%
%   USAGE :
%       [R, dRdU, m] = STATRESFUN(m, U, Fstat)
%   INPUTS :
%       U       : (Nd, 1) vector of DOFS
%       Fstat   : (Nd, 1) vector of inhomogeneous forcing
%   OUTPUTS :
%       R       : (Nd, 1) Residue vector (KU+FNL-Fstat)
%       dRdU    : (Nd, Nd) Jacobian matrix (K+JNL)

%     [FNL, dFNL, ~, m] = m.NLFORCE(0, U, U*0, -1, 1);
    [FNL, dFNL, ~, m] = m.NLFORCE(0, U, U*0, 0, 1);
    
    % Residue
    R = m.K*U+FNL-Fstat;
    dRdU = m.K+dFNL;
end