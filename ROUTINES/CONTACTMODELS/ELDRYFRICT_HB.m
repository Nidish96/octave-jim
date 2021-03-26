function [FXYN, JXYN] = ELDRYFRICT_HB(UXYN, h, Nt, kxynmu, N0, varargin)
%ELDRYFICT returns the forces and jacobians for the elastic dry friction model in the frequency domain
%
% USAGE:
% ------  
%   [FXYN, JXYN] = ELDRYFRICT(UXYN, h, Nt, kxynmu, N0);
% INPUTS:
% -------
%   UXYN	  : (Np*3*Nhc,1)
%	 [[ux1; uy1; uz1; ux2; ...]_A0; [ux1; uy1; uz1; ux2; ...]_A1; [ux1; uy1; uz1; ux2; ...]_B1; ...]
%   h		  : List of harmonics
%   Nt		  : Number of time points for AFT
%   kxynmu	  : (Np*4,1) [[kx1; ky1; kn1; mu1]; [kx2; ky2; kn2; mu2]; ...]
%   N0		  : (1,Np) [N0]
% OUTPUTS:
% --------
%   FXYN	  : (Np*3*Nhc,1)
%	 [[fx1; fy1; fz1; fx2; ...]_A0; [fx1; fy1; fz1; fx2; ...]_A1; [fx1; fy1; fz1; fx2; ...]_B1; ...]
%   JXYN	  : (Np*3*Nhc, Np*3*Nhc)  

    error('This is work in progress and has to be debugged');

  if length(varargin)>=1
    tol = varargin{1};
  else
    tol = single(1e-6);
  end
  
  Nhc = uint32(sum(h==0)+2*sum(h~=0));
  Np = uint32(length(UXYN)/(3*Nhc));

  fxyn = zeros(Np*3, Nt, 'single');
  jxyn = zeros(Nt, Np*3, Np*3*Nhc, 'single');

  cst = TIMESERIES_DERIV(Nt, h, eye(Nhc, 'single'), 0);  % Fourier-Galerkin Bases

  uxyn = TIMESERIES_DERIV(Nt, h, reshape(UXYN, Np*3, Nhc)', 0)';  % 3Np x Nt

  % TIME DOMAIN CALCULATIONS
  % NORMAL FORCE AND JACOBIAN
  fxyn(3:3:end, :) = max(repmat(kxynmu(3:4:end),1,Nt).*uxyn(3:3:end,:)+N0, 0);
  for p=uint32(1:Np)
    ip0 = (p-1)*3;
    for t=uint32(1:Nt)
      if fxyn(ip0+3, t)~=0  % point in contact
	jxyn(t, ip0+3, (ip0+3):(3*Np):(3*Np*Nhc)) = kxynmu((p-1)*4+3)*cst(t, :);
      end
    end
  end

  convstat = repmat(false, Np, 1);  % Convergence status
  while any(~convstat) % while some point hasn't converged
    fprev = fxyn(:, 1);  % First point in time from previous estimates

    for p=uint32(1:Np)
				% Check if point already converged
      if convstat(p)
	continue;
      end
			       % Conduct calculations if not converged
      ip0 = double((p-1)*3);
      for t=uint32(1:Nt)
	tm1 = mod(t-2,Nt)+1;
				% CHECK FOR SEPARATION
	if fxyn(ip0+3, t)==0
	  fxyn(ip0+(1:2), t) = 0;
	  jxyn(t, ip0+(1:2), ip0+1:(3*Np):(3*Np*Nhc)) = 0;
	  jxyn(t, ip0+(1:2), ip0+2:(3*Np):(3*Np*Nhc)) = 0;
	  jxyn(t, ip0+(1:2), ip0+3:(3*Np):(3*Np*Nhc)) = 0;
	  continue;
	end
				% STICK PREDICTION
	fxyn(ip0+1, t) = kxynmu((p-1)*4+1)*(uxyn(ip0+1,t)-uxyn(ip0+1,tm1)) + ...
			 fxyn(ip0+1, tm1);  % fx
	fxyn(ip0+2, t) = kxynmu((p-1)*4+2)*(uxyn(ip0+2,t)-uxyn(ip0+2,tm1)) + ...
			 fxyn(ip0+2, tm1);  % fy
	
	jxxstick = permute(kxynmu((p-1)*4+1)*(cst(t,:)-cst(tm1,:)), [1, 3, 2]) + ...
		   jxyn(tm1, ip0+1, ip0+1:(3*Np):(3*Np*Nhc));
	jxystick = jxyn(tm1,ip0+1, ip0+2:(3*Np):(3*Np*Nhc));
	jxnstick = jxyn(tm1,ip0+1, ip0+3:(3*Np):(3*Np*Nhc));

	jyxstick = jxyn(tm1,ip0+2, ip0+1:(3*Np):(3*Np*Nhc));
	jyystick = permute(kxynmu((p-1)*4+2)*(cst(t,:)-cst(tm1,:)), [1, 3, 2]) + ...
		   jxyn(tm1,ip0+2, ip0+2:(3*Np):(3*Np*Nhc));
	jynstick = jxyn(tm1,ip0+2, ip0+3:(3*Np):(3*Np*Nhc));
	
				% CHECK FOR SLIP
	fT = sqrt(fxyn(ip0+1, t)^2+fxyn(ip0+2, t)^2);
	fslip = kxynmu((p-1)*4+4)*fxyn(ip0+3, t);  % mu*fN
	if fT<fslip  % Point is stuck
		     % jxx
	  jxyn(t, ip0+1, ip0+1:(3*Np):(3*Np*Nhc)) = jxxstick;
				% jxy
	  jxyn(t, ip0+1, ip0+2:(3*Np):(3*Np*Nhc)) = jxystick;
				% jxn
	  jxyn(t, ip0+1, ip0+3:(3*Np):(3*Np*Nhc)) = jxnstick;

				% jxz
	  jxyn(t, ip0+2, ip0+1:(3*Np):(3*Np*Nhc)) = jyxstick;
				% jyy
	  jxyn(t, ip0+2, ip0+2:(3*Np):(3*Np*Nhc)) = jyystick;
				% jyz
	  jxyn(t, ip0+2, ip0+3:(3*Np):(3*Np*Nhc)) = jynstick;
	else % Point has slipped
	  jxyn(t, ip0+1, ip0+1:(3*Np):(3*Np*Nhc)) = fslip*(fxyn(ip0+2, t)^2*jxxstick-prod(fxyn(ip0+(1:2), t))*jyxstick)/fT^(3); %jxx
	  jxyn(t, ip0+1, ip0+2:(3*Np):(3*Np*Nhc)) = fslip*(fxyn(ip0+2, t)^2*jxystick-prod(fxyn(ip0+(1:2), t))*jyystick)/fT^(3); %jxy
	  jxyn(t, ip0+1, ip0+3:(3*Np):(3*Np*Nhc)) = kxynmu((p-1)*4+4)*jxyn(t, ip0+3, ip0+3:(3*Np):(3*Np*Nhc))*fxyn(ip0+1, t)/fT + ...
						    fslip*(fxyn(ip0+2, t)^2*jxnstick-prod(fxyn(ip0+(1:2), t))*jynstick)/fT^(3); %jxn

	  jxyn(t, ip0+2, ip0+1:(3*Np):(3*Np*Nhc)) = fslip*(fxyn(ip0+1, t)^2*jyxstick-prod(fxyn(ip0+(1:2), t))*jxxstick)/fT^(3); %jyx
	  jxyn(t, ip0+2, ip0+2:(3*Np):(3*Np*Nhc)) = fslip*(fxyn(ip0+1, t)^2*jyystick-prod(fxyn(ip0+(1:2), t))*jxystick)/fT^(3); %jyy
	  jxyn(t, ip0+2, ip0+3:(3*Np):(3*Np*Nhc)) = kxynmu((p-1)*4+4)*jxyn(t, ip0+3, ip0+3:(3*Np):(3*Np*Nhc))*fxyn(ip0+2, t)/fT + ...
						    fslip*(fxyn(ip0+1, t)^2*jynstick-prod(fxyn(ip0+(1:2), t))*jxnstick)/fT^(3); %jyn

	  fxyn(ip0+1, t) = fslip*fxyn(ip0+1, t)/fT; % fx
	  fxyn(ip0+2, t) = fslip*fxyn(ip0+2, t)/fT; % fy
	end
      end % Iterations over quadrature points
    end  % Iterations over time

    tmp = fxyn(:, 1);  tmp(tmp==0) = 1;  % relative scaling for checking convergence
    convstat = sum(reshape(abs(fxyn(:, 1) - fprev)./tmp, 3, Np), 1)<tol;
  end  % All points have converged
				% FREQUENCY DOMAIN TRANSFORMATION
  FXYN = reshape(GETFOURIERCOEFF(h, fxyn')', 3*Np*Nhc, 1);
  clear fxyn
  
  % % % THIS CONSUMES A LOT OF RAM - BETTER TO IMPLEMENT THE FFT IN-PLACE - AND DON'T USE PERMUTE
  % jxyn = permute(jxyn, [1 3 2]);
  % JXYN = GETFOURIERCOEFF(h, reshape(jxyn, Nt, 3*Np*Nhc*3*Np));
  % JXYN = reshape(JXYN', Np*3*Nhc, Np*3*Nhc)';

				% USING LOOPS INSTEAD OF PERMUTE
  jxyn = GETFOURIERCOEFF(h, reshape(jxyn, Nt, 3*Np*Nhc*3*Np));
  jxyn = reshape(jxyn, Nhc, Np*3, Np*3*Nhc);
  JXYN = zeros(Np*3*Nhc, Np*3*Nhc, 'single');
  for ih=1:Nhc
    for jh=1:Nhc
      JXYN((ih-1)*Np*3+(1:Np*3), (jh-1)*Np*3+(1:Np*3)) = ...
      reshape(jxyn(ih, :, (jh-1)*Np*3+(1:Np*3)), Np*3, Np*3);
    end
  end
end
