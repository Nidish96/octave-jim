function [U, R, eflag, it, Jc] = NSOLVE(func, U0, opts)
%NSOLVE Uses Newton iterations to solve
%
% USAGE:
%  [U, R, eflag, it, jc] = NSOLVE(func, U0, opts);
% INPUTS:
%  func		: Function handle [F, dFdU] = func(u);
%  U0		: Initial Guess
%  opts		: Options structure with,
% 	reletol (float)
% 	etol    (float)
% 	utol    (float)
% 	rtol    (float)
% 	Display (boolean)
% OUTPUTS:
%  U		:
%  R		:
%  eflag	:
%  it		:
%  Jc		:

  if ~isfield(opts, 'Dscale')
    opts.Dscale = ones(size(U0));
  end
  Nu = length(U0);
  
  [R0, J0] = func(opts.Dscale.*U0);
  dU0 = (-J0\R0)./opts.Dscale;
  e0  = abs(R0'*dU0);
  
  if (e0 < eps)
    e0 = 1.0;
    R0 = ones(size(R0));
  end
  r0 = sqrt(mean(R0.^2));
  u0 = sqrt(mean(dU0.^2));
  
  r   = r0;
  e   = e0;
  u   = u0;

  U  = U0;
  dU = dU0;
  it = 0;

  eflag = 8*(e<opts.etol) + 4*(e/e0<opts.reletol) + 2*(r<opts.rtol) + 1*(u<opts.utol) + 16*(e/eps<1e3);
  if opts.Display
    fprintf('ITN, E, E/E0, r, du\n%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, eflag);
  end
  eflagp = 0;
  
  while (eflag<6 || it==0) && it<=opts.ITMAX
    U  = U + dU;
    it = it+1;
    
    [R, Jc] = func(opts.Dscale.*U);
    dU = (-Jc\R)./opts.Dscale;
    
    e = abs(R'*dU);
    r = sqrt(mean(R.^2));
    u = sqrt(mean(dU.^2));

    eflag = 8*(e<opts.etol) + 4*(e/e0<opts.reletol) + 2*(r<opts.rtol) + 1*(u<opts.utol) + 16*(e/eps<1e3);
    if opts.Display
      fprintf('%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, eflag);
    end
    
    if it>opts.ITMAX
      eflag = 0;
      break;
    end
  end

  % Rescale Solution
  U = opts.Dscale.*U;
  
  if eflag == 0
    disp('No Convergence : Returning')
    return;
  end

  if opts.Display
    disp('Successful Campaign!')
  end
end
