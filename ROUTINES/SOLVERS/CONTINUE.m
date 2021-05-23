function [U, dUdlam, Ss] = CONTINUE(func, u0, lam0, lam1, ds, varargin)
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
% error('This is obsolete now - please use PRECOCONT.m')
  %% Default options
  Copt = struct('Nmax', 100, 'dsmax', ds*5, 'dsmin', ds/5, ...
                'angopt', pi/6, 'startdir', sign(lam1-lam0),...
                'Display', true, 'itDisplay', false, 'lsrch', 0, ...
                'nev', 1, 'adj', 1, 'itopt', 10, ...
		'arclengthparm', 'orthogonal', 'ITMAX', 20, ...
        'Dscale', ones(length(u0)+1,1), 'DynDscale', 0, ...
                'opts', struct('reletol', 1e-6, 'rtol', 1e-6, 'etol', ...
                              1e-6, 'utol', 1e-6, ...
                              'Display', false, 'lsrch', 0));
  if nargin==6
    nflds = fieldnames(varargin{1});
    for i=1:length(nflds)
      Copt.(nflds{i}) = varargin{1}.(nflds{i});
    end
  end
  Copt.opts.lsrch = Copt.lsrch;
  Copt.opts.Display = Copt.itDisplay;
  Copt.opts.ITMAX = Copt.ITMAX;

  %% Allocations
  U = zeros(length(u0)+1, Copt.Nmax);
  dUdlam = zeros(length(u0)+1, Copt.Nmax);

  J0 = zeros(length(u0)+1);  
  
  ds0 = ds;
  
  %% Correct initial solution
  Copt.opts.Dscale = Copt.Dscale(1:end-1);
  tm = Copt.opts.Display;
  Copt.opts.Display = true;
  [U(1:end-1,1), ~, eflag] = NSOLVE(@(u) func([u; lam0]), u0, Copt.opts);
  Copt.opts.Display = tm;
  U(end, 1) = lam0;
  if eflag < 0
      error('Initial point non-convergent!');
  elseif Copt.Display
      disp('Initial Point Converged')
      fprintf('Starting Continuation from %f to %f\n', lam0, lam1);
  end
  
  if Copt.DynDscale
      Copt.opts.Dscale = max(abs(U(:, 1)), Copt.Dscale);
  else
      Copt.opts.Dscale = Copt.Dscale;
  end
  
  arclenparms = {'arclength', 'orthogonal'};
  
  % Extract gradient information
  [~, J0(1:end-1,1:end-1), J0(1:end-1,end)] = func(U(:, 1));
  J0 = J0*diag(Copt.opts.Dscale);
  dUdlam(:, 1) = [-J0(1:end-1,1:end-1)\J0(1:end-1,end); 1];
  
  % Initial tangent
  dxn = Copt.startdir;
  z = dUdlam(1:end-1, 1);
  al = dxn/sqrt(1+z'*z);
  alp = al;
  Ss = zeros(Copt.Nmax,1);
  
  %% BEGIN CONTINUATION
  lam  = lam0;  % current parameter
  lamp = lam0;  % previous parameter
  n    = 1;
  
  tries = 0;
  
%   oopts = optimset('Jacobian', 'on', 'Display', 'off');  % options for fsolve if NSOLVE fails 
  
  dUdlam(:, 1) = Copt.opts.Dscale.*dUdlam(:, 1);
  
  dUn = [z;1];
  u0 = U(:, 1) + ds*(al*dUn);
  while ( (lam-lam1)*(lamp-lam1) >= 0 && n<Copt.Nmax )
      if Copt.DynDscale
          Copt.opts.Dscale = max(abs(U(:, n+1)), Copt.opts.Dscale);
      end
        
%     [U(:, n+1), ~, eflag, op, J0] = fsolve(@(u) EXRES(func, u, U(:, n), al*dUdlam(:, n), ds, Copt.arclengthparm), u0, oopts);
%     its = op.iterations;      
%     [U(:, n+1), ~, eflag, its, J0] = NSOLVE(@(u) EXRES(func, u, U(:, n), al*dUdlam(:, n), ds, Copt.arclengthparm), u0, Copt.opts);

%     [U(1:end-1, n+1), ~, eflag, op] = fsolve(@(u) ELRES(func, u, U(:, n), al*dUdlam(:, n), ds, Copt.arclengthparm), u0(1:end-1), oopts);
%     its = op.iterations;      
%     [U(1:end-1, n+1), ~, eflag, its] = NSOLVE(@(u) ELRES(func, u, U(:, n), al*dUdlam(:, n), ds, Copt.arclengthparm), u0(1:end-1), Copt.opts);

    [U(:, n+1), ~, eflag, its, J0, dUn] = ELIMSOLVE(func, u0, U(:, n), al*dUn, ds, Copt.arclengthparm, Copt.opts);
    dUn = dUn/al;
    if eflag<=0
        fprintf('Diverged %d - ', tries);
        switch tries
            case {0}
                if ds>Copt.dsmin
                    fprintf('reducing step-size\n');
                    ds = ds/2;
                    tries = tries-1;
                elseif ds~=ds0
                    fprintf('resetting step-size\n')
                    ds = ds0;
                end
                u0 = U(:, n) + ds*al*dUn;
            case {2, 3}
                fprintf('trying eigenvector update\n');
                
%                 [~, dR0, lnx, J0] = ELRES(func, u0(1:end-1), U(:, n), al*dUn, ds, Copt.arclengthparm);
                [eV, eD] = eigs(dR0, 10, 'SM');
                [eD, si] = sort(diag(eD));
                eV = eV(:, si); 
                
                eU = -eV(:, 1)/norm(eV(:, 1));
                
                if tries==2
                    u0 = [U(1:end-1, n) + eU*ds; lnx];
                else
                    u0 = [U(1:end-1, n) + -eU*ds; lnx];
                end
            case {3}
                fprintf('retrying with Branch starter\n')
                u0 = U(:, 1);
            case {4}
                fprintf('retrying with zeros\n')
                u0 = U(:, 1)*0;
            case {5}
                fprintf('maxing step-size\n')
                ds = Copt.dsmax;
                u0 = U(:, n) + ds*al*dUn;
            otherwise
                fprintf('Stopping\n');
                break;
        end
        tries = tries+1;
        continue;
    else
%         [~, dR0, U(end, n+1), J0] = ELRES(func, U(1:end-1, n+1), U(:, n), al*dUn, ds, Copt.arclengthparm);
        tries = 0;
    end
    
    Ss(n+1) = Ss(n)+ds;
    J0 = J0*diag(Copt.opts.Dscale);
    dUdlam(:, n+1) = [-J0(1:end-1,1:end-1)\J0(1:end-1,end); 1];
    dxn = sign(dxn*dUn'*dUdlam(:, n+1));
%     dxn = sign((U(end,n+1)-U(end,n))*dUdlam(:, n)'*dUdlam(:, n+1));

				% Step Update
    lamp = lam;
    lam  = U(end, n+1);
				% tangent and predictor
    z = dUdlam(1:end-1, n+1);
    al = dxn/sqrt(1+z'*z);
    
    dUdlam(:, n+1) = Copt.opts.Dscale.*dUdlam(:, n+1);
    
				% step size adaptation
    theta = acos(alp*al*dUn'*[z; 1]);
    if Copt.Display
      fprintf('%d %f %f %e %d (%d)\n', n+1, U(end,n+1), ds, theta, its, eflag);
    end
    if n>1
        if abs(theta)>1.1*Copt.angopt % angle check
            ds = max(Copt.dsmin, ds/2);
        elseif abs(theta)<0.9*Copt.angopt  % angle and convergence check
            ds = min(max(ds*sqrt(Copt.itopt/its), Copt.dsmin), Copt.dsmax);
        end
%         ds = min(max(ds*sqrt(Copt.itopt/its), Copt.dsmin), Copt.dsmax);
    end
    alp = al;    
    n = n+1;
    
    dUn = [z; 1];

    u0 = U(:, n) + ds*al*dUn;
  end
  U = U(:, 1:n);
  dUdlam = dUdlam(:, 1:n);
  Ss = Ss(1:n);
  
  if (lam-lam1)*(lamp-lam1) < 0
    disp('Continuation completed successfully');
  else 
    disp('Premature Termination');
  end
end

%%
function [ul, Re, eflag, its_tot, Je, ulp0] = ELIMSOLVE(func, ul, ul0, ulp0, ds, parm, opts)
%ELIMSOLVE is the residual function along with the continuation constraint
%(with elimination)
%  USAGE:
%   [Re, Je] = ELIMSOLVE(func, ul, ul0, ulp0, ds, parm, opts)
%  INPUTS:
%   ul   : Solution vector (with param) - (Nd, 1)
%   ul0 : previous point (with param) - (Nd+1, 1)
%   ulp0: previous tgt (with param) - (Nd+1, 1)
%   ds  : step
%   parm: parameterization
%  OUTPUTS:
%   Re  : (Nd, 1)
%   Je  : (Nd, Nd)
%   dRe : (Nd+1, Nd+1)

  l0 = ul0(end);
  lp = ul(end);

  sgn = sign(ulp0(end));

  swch = 0;
  
  its_tot = 0;
  while true
      switch parm
          case {'orthogonal' , 'Orthogonal'}
              l = l0 + (ds-ulp0(1:end-1)'*(ul(1:end-1)-ul0(1:end-1)))/ulp0(end);
              dcR = ulp0';
          case {'arclength' , 'Arclength'}
              sgn = sign(ulp0(end));              
              l = l0 + sgn*sqrt(ds^2 - (ul(1:end-1)-ul0(1:end-1))'*(ul(1:end-1)-ul0(1:end-1)));

              if ~isreal(l)
                  l = lp;
              else
                dcR = 2*(ul-ul0)';
              end
          otherwise
              sgn = sign(ulp0(end));
              [Re0, dRdUe0, dRdLe0] = func(ul0);

              b = 0.5;
              c = (1-b)/(dRdLe0'*dRdUe0*dRdLe0);

              du = (ul-ul0(1:end-1));
              l = l0 + sgn*sqrt(ds^2 - c*(du'*dRdUe0*du))/b;
              
              dcR = [2*c*(du'*dRdUe0) 2*b*(l-l0)];
      end

      % Corrector 
      [un, ~, eflag, its] = NSOLVE(@(u) func([u; l]), ul(1:end-1), opts);
      
      if eflag<=0  % Try direct solve
        [unl, ~, eflag, its] = NSOLVE(@(ul) EXRES(func, ul, ul0, ulp0, ds, parm), ul0+ulp0*ds/2, opts);
        if eflag<=0
            if swch==1
                break;
            else
                swch = 1;
                ulp0 = -ulp0;
                continue;
            end
        end
        swch = 0;
            
        un = unl(1:end-1);
        l = unl(end);
        lp = l;  % no repeats necessary if direct solve works
        sgn = sign(l-l0);
      end

      its_tot = its_tot + its;
%       its_tot = its_tot + 1;
      
      % Update
      ul = [un; l];
      if abs((l-lp)/lp)<opts.reletol
          break;
      else
          lp = l;
      end
  end
  
  [Re, dRdUe, dRdLe] = func(ul);
  Je = [dRdUe dRdLe; dcR];
end

%%
function [Re, Je, l, dRe] = ELRES(func, u, ul0, ulp0, ds, parm)
%EXRES is the residual function along with the continuation constraint
%(with elimination)
%  USAGE:
%   [Re, Je] = ELRES(func, u, ul0, ulp0, ds, parm)
%  INPUTS:
%   u   : Solution vector (no param) - (Nd, 1)
%   ul0 : previous point (with param) - (Nd+1, 1)
%   ulp0: previous tgt (with param) - (Nd+1, 1)
%   ds  : step
%   parm: parameterization
%  OUTPUTS:
%   Re  : (Nd, 1)
%   Je  : (Nd, Nd)
%   dRe : (Nd+1, Nd+1)

  l0 = ul0(end);
  [Re0, dRdUe0, dRdLe0] = func(ul0);
  switch parm
      case {'orthogonal' , 'Orthogonal'}
          l = l0 + (ds-ulp0(1:end-1)'*(u-ul0(1:end-1)))/ulp0(end);
          dldu = -ulp0(1:end-1);
      case {'arclength' , 'Arclength'}
          sgn = sign(ulp0(end));
          l = l0 + sgn*sqrt(ds^2 - (u-ul0(1:end-1))'*(u-ul0(1:end-1)));
          
          if ~isreal(l)
              disp('hey');
          end
          
          dldu = -(u-ul0(1:end-1))'/(l-l0);
      otherwise
          sgn = sign(ulp0(end));
          b = 0.5;
          c = (1-b)/(dRdLe0'*dRdUe0*dRdLe0);
    
          du = (u-ul0(1:end-1));
          l = l0 + sgn*sqrt(ds^2 - c*(du'*dRdUe0*du))/b;
          
          dldu = -c*du'*dRdUe0/(b*(l-l0));
  end

  [Re, dRdUe, dRdLe] = func([u; l]);
  Je = dRdUe + dRdLe*dldu;
  
  if nargout>=4
      switch parm
        case {'orthogonal' , 'Orthogonal'}
    %       Re = [Re; up0'*(u-u0)-ds];
          dRe = [dRdUe dRdLe; ulp0'];
        case {'arclength' , 'Arclength'}
    %       Re = [Re; (u-u0)'*(u-u0)-ds^2];
          dRe = [dRdUe dRdLe; 2*([u; l]-ul0)'];
      end
  end
end

%%
function [Re, Je] = EXRES(func, u, u0, up0, ds, parm)
%EXRES is the residual function along with the continuation constraint
%(without elimination)

  [Re0, dRdUe0, dRdLe0] = func(u0);
  [Re, dRdUe, dRdLe] = func(u);
  switch parm
    case {'orthogonal' , 'Orthogonal'}
      Re = [Re; up0'*(u-u0)-ds];
      Je = [dRdUe dRdLe; up0'];
    case {'arclength' , 'Arclength'}
      Re = [Re; (u-u0)'*(u-u0)-ds^2];
      Je = [dRdUe dRdLe; 2*(u-u0)'];
    otherwise        
      b = 0.5;
      c = (1-b)/(dRdLe0'*dRdUe0*dRdLe0);
    
      du = (u(1:end-1)-u0(1:end-1));
      dL = (u(end)-u0(end));
      Re = [Re; 
          c*(du'*dRdUe0*du)+b*dL^2-ds^2];
      Je = [dRdUe dRdLe;
          c*(2*du'*dRdUe0), 2*b*dL];
  end
end
