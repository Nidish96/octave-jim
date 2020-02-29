function [R, dRdU, dRda, Z, dRdPar] = QS_RESFUN(Ua, Z, Pars, L, pA, MESH, M, K, Fs, Fv)
%QS_RESFUN returns the residual function for RQNMA
%
%  USAGE:
%  ------  
%    [R, Z, dRdU, dRda, dRdPar] = QS_RESFUN(Ua, Z, Pars, L, pA, MESH, M, K, Fs, Fv);
%  INPUTS:
%  -------
%    Ua		: (Nu+1,1) [U; a]
%    Z		: (Nh1,Nh2)
%    Pars	: (Npars,1)
%    L		: (Nph,Nu) Projection/Mask matrix for U
%    pA		: (Nparsp,Npars) Projection/Mask matrix for Pars
%    MESH	: MESH2D class object
%    M		: (Nu,Nu) Mass Matrix
%    K		: (Nu,Nu) Stiffness Matrix
%    Fs		: (Nu,1) Static Forcing Amplitude
%    Fv		: (Nu,1) Variable Forcing Amplitude (amplified by a)
  
  % Compute Contact Forces
  [Fnl, Z, dFdUnl, ~, dFnldP] = CONTACTEVAL(MESH, Ua(1:end-1), Z, Ua(1:end-1)*0, Pars, pA, L);

  % Compute Residuals
  R = K*Ua(1:end-1) + Fnl - Ua(end)*Fv - Fs;
  dRdU = K + dFdUnl;
  dRda = -Fv;
  dRdPar = dFnldP;
end
