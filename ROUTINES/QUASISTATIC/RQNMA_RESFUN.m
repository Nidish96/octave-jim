function [R, dRdUl, dRdq, Z, dRdPar] = RQNMA_RESFUN(Ulq, Z, Pars, L, pA, MESH, M, K, Fs, Us)
%RQNMA_RESFUN returns the residual function for RQNMA
%
%  USAGE:
%  ------  
%    [R, Z, dRdUl, dRdq, dRdPar] = RQNMA_RESFUN(Ulq, Z, Pars, L, pA, MESH, M, K, Fs, Us);
%  INPUTS:
%  -------
%    Ulq	: (Nu+2,1)
%    Z		: (Nh1,Nh2)
%    MESH	: MESH2D class object
%    M		: (Nu,Nu) Mass Matrix
%    K		: (Nu,Nu) Stiffness Matrix
%    L		: (Nph,Nu) Projection/Mask matrix
%    Fs		: (Nu,1) Static Forcing Amplitude
%    Us		: (Nu,1) Guess for Static solution

  % Compute Contact Forces
  [Fnl, Z, dFdUnl, ~, dFnldP] = CONTACTEVAL(MESH, Ulq(1:end-2), Z, Ulq(1:end-2)*0, Pars, pA, L);

  % % Compute Residuals
  % R = [K*Ulq(1:end-2)+Fnl-Fs-Ulq(end-1)*M*(Ulq(1:end-2)-Us);
  %      0.5*((Ulq(1:end-2)-Us)'*M*(Ulq(1:end-2)-Us)-Ulq(end)^2)];
  % dRdUl = [K+dFdUnl-Ulq(end-1)*M, -M*(Ulq(1:end-2)-Us);
  % 	   (Ulq(1:end-2)-Us)'*M, 0];
  % dRdq = [zeros(size(Ulq(1:end-2))); Ulq(end)];
  % dRdPar = [dFnldP; zeros(1, length(Pars))];

  % % Compute Residuals - Improve Conditioning by multiplying second equation by lambda?
  % R = [K*Ulq(1:end-2)+Fnl-Fs-Ulq(end-1)*M*(Ulq(1:end-2)-Us);
  %      Ulq(end-1)*0.5*((Ulq(1:end-2)-Us)'*M*(Ulq(1:end-2)-Us)-Ulq(end)^2)];
  % dRdUl = [K+dFdUnl-Ulq(end-1)*M, -M*(Ulq(1:end-2)-Us);
  % 	   Ulq(end-1)*(Ulq(1:end-2)-Us)'*M, 0.5*((Ulq(1:end-2)-Us)'*M*(Ulq(1:end-2)-Us)-Ulq(end)^2)];
  % dRdq = [zeros(size(Ulq(1:end-2))); Ulq(end-1)*Ulq(end)];
  % dRdPar = [dFnldP; zeros(1, length(Pars))];

  % % Compute Residuals - Improve Conditioning by rewriting in terms of omega (sqrt(lambda))
  % R = [K*Ulq(1:end-2)+Fnl-Fs-Ulq(end-1)^2*M*(Ulq(1:end-2)-Us);
  %      Ulq(end-1)*0.5*((Ulq(1:end-2)-Us)'*M*(Ulq(1:end-2)-Us)-Ulq(end)^2)];
  % dRdUl = [K+dFdUnl-Ulq(end-1)^2*M, -2*Ulq(end-1)*M*(Ulq(1:end-2)-Us);
  % 	   Ulq(end-1)*(Ulq(1:end-2)-Us)'*M, 0.5*((Ulq(1:end-2)-Us)'*M*(Ulq(1:end-2)-Us)-Ulq(end)^2)];
  % dRdq = [zeros(size(Ulq(1:end-2))); Ulq(end-1)*Ulq(end)];
  % dRdPar = [dFnldP; zeros(1, length(Pars))];

  % Compute Residuals - Combine equations and substitute constraint once
  % SLIGHTLY BETTER CONDITIONING
  R = [K*Ulq(1:end-2)+Fnl-Fs-Ulq(end-1)*M*(Ulq(1:end-2)-Us);
       (Ulq(1:end-2)-Us)'*(K*Ulq(1:end-2)+Fnl-Fs)-Ulq(end-1)*Ulq(end)^2];

  dRdUl = [K+dFdUnl-Ulq(end-1)*M, -M*(Ulq(1:end-2)-Us);
	   Ulq(1:end-2)'*K+Fnl'-Fs'+(Ulq(1:end-2)-Us)'*(K+dFdUnl), -Ulq(end)^2];
  dRdq = [zeros(size(Ulq(1:end-2))); -2*Ulq(end-1)*Ulq(end)];
  dRdPar = [dFnldP; (Ulq(1:end-2)-Us)'*dFnldP];

  % NEED TO COME UP WITH ELIMINATED CONSTRAINT FORMULATION TO BE BETTER OFF
end
