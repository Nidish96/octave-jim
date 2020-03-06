function [U, dUdlam] = CONTINUE(func, u0, lam0, lam1, ds, varargin)
%CONTINUE Conducts the continuation and solves the system
%
% USAGE:
%   [U, dUdlam] = CONTINUE(func, u0, lam0, lam1, ds, varargin);
% INPUTS:
%   func	: Function handle (Jacobian is strictly necessary)
%   u0		: (Nu,1) initial guess for initial point
%   lam0	: Parameter value at initial point
%   lam1	: Parameter value at final point
%   ds		: Arc length parameter at initial point
%   Copt 	: (optional) Options structure
% OUTPUTS:
%   U		: (Nu+1, Np) Solution vector
%   dUdlam	: (Nu+1, Np) Solution derivative wrt parameter

				% Default options
  Copt = struct('Nmax', 100, 'dsmax', ds*5, 'dsmin', ds/5, ...
                'angopt', pi/6, 'startdir', sign(lam1-lam0),...
                'Display', true, 'nev', 1, 'adj', 1,...
		'arclengthparm', 'orthogonal', ...
                'opts', struct('reletol', 1e-6, 'rtol', 1e-6, 'etol', ...
                              1e-6, 'utol', 1e-6, 'ITMAX', 20, ...
                              'Display', true));
  if nargin==6
    nflds = fieldnames(varargin{1});
    for i=1:length(nflds)
      Copt.(nflds{i}) = varargin{1}.(nflds{i});
    end
  end

  % Allocations
  U = zeros(length(u0)+1, Copt.Nmax);
  dUdlam = zeros(length(u0)+1, Copt.Nmax);

  R0 = zeros(length(u0)+1, 1);
  J0 = zeros(length(u0)+1);  
  
  % Correct initial solution
  if isfield(Copt, 'Dscale')
      Copt.opts.Dscale = Copt.Dscale(1:end-1);
  end
  [U(1:end-1,1), ~, eflag] = NSOLVE(@(u) func([u; lam0]), u0, Copt.opts);
  U(end, 1) = lam0;
  if eflag < 0
      error('Initial point non-convergent!');
  elseif Copt.Display
      disp('Initial Point Converged')
  end
  % Extract gradient information
  [~, J0(1:end-1,1:end-1), J0(1:end-1,end)] = func(U(:, 1));
  dUdlam(:, 1) = [-J0(1:end-1,1:end-1)\J0(1:end-1,end); 1];
  
				% BEGIN CONTINUATION
  lam  = lam0;  % current parameter
  lamp = lam0;  % previous parameter
  n    = 1;
  
				% Initial tangent
  dxn = Copt.startdir;
  al = dxn/sqrt(1+sum(dUdlam(1:end-1, 1).^2));
  
  alp = al;
  
  u0 = U(:, 1) + ds*al*dUdlam(:, 1);
  if isfield(Copt, 'Dscale')
      Copt.opts.Dscale = Copt.Dscale;
      clear Copt.Dscale
  end
  while ( (lam-lam1)*(lamp-lam1) >= 0 && n<Copt.Nmax )
    [U(:, n+1), R0, eflag, its, J0] = NSOLVE(@(u) EXRES(func, u, U(:, n), al*dUdlam(:, n), ds, Copt.arclengthparm), u0, Copt.opts);
    if eflag<=0
      if ds <= Copt.dsmin
	disp('Diverged! Trying with first initial guess!');
	[U(:, n+1), R0, eflag, its, J0] = NSOLVE(@(u) EXRES(func, u, U(:, n), al*dUdlam(:, n), ds, Copt.arclengthparm), U(:, 1), Copt.opts);
	if eflag<=0
	  disp('Diverged! Returning');
	  break;
	end
      else
        disp('Diverged - reducing step-size');
        ds = ds/2;
        u0 = U(:, n) + ds*duds;
        continue;
      end
    end
    dUdlam(:, n+1) = [-J0(1:end-1,1:end-1)\J0(1:end-1,end); 1];
    dxn = sign(dxn*sum(prod(dUdlam(:, n:n+1),2)));
				% Step Update
    lamp = lam;
    lam  = U(end, n+1);
				% tangent and predictor
    al = dxn/sqrt(1+sum(dUdlam(1:end-1, n+1).^2));   
    
				% step size adaptation
    theta = acos(alp*al*sum(prod(dUdlam(:, n:n+1),2)));
    if abs(theta)>1.1*Copt.angopt % angle check
      ds = max(Copt.dsmin, ds/2);
    elseif its<=10 && abs(theta)<0.9*Copt.angopt  % angle and convergence check
      ds = min(Copt.dsmax, ds*2);
    end
    if Copt.Display
      fprintf('%d %f %f %e %d\n', n+1, U(end,n+1), ds, theta, eflag);
    end
    
    alp = al;    
    n = n+1;

    u0 = U(:, n) + ds*al*dUdlam(:, n);
  end
  U = U(:, 1:n);
  dUdlam = dUdlam(:, 1:n);
  
  if (lam-lam1)*(lamp-lam1) < 0
    disp('Continuation completed successfully');
  else 
    disp('Premature Termination');
  end
end

function [Re, Je] = EXRES(func, u, u0, up0, ds, parm)
%EXRES is the residual function along with the continuation constraint
  [Re, dRdUe, dRdLe] = func(u);
  switch parm
    case {'orthogonal' , 'Orthogonal'}
      Re = [Re; up0'*(u-u0)-ds];
      Je = [dRdUe dRdLe; up0'];
    case {'arclength' , 'Arclength'}
      Re = [Re; (u-u0)'*(u-u0)-ds^2];
      Je = [dRdUe dRdLe; 2*(u-u0)'];
  end
end