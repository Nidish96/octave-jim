function [fxyn, DfxynDuxyn] = ELDRYFRICT(t, uxyn, Kt, kn, mu, gap, varargin)
%ELDRYFICT returns the forces and jacobians for the elastic dry friction model
%
% USAGE:
% ------  
%   [fxyn, zxy, DfxynDuxyn] = ...
%      ELDRYFRICT(t, uxyn, Kt, kn, mu, gap, tp, uxynp, fxynp, dfxynp, h)
% INPUTS:
% -------
%   uxyn	  : (3*Np, 1) [[ux; uy; un]_1; [ux; uy; un]_2; ...]
%   Kt		  : (3, Np) [[kx; ky; kxy]_1, [kx; ky; kxy]_2, ...]
%   kn 		  : (1, Np) [kn1, kn2, ...]
%   mu 		  : (1, Np) [mu1, ...]
%   gap	 	  : (1, Np)
% OUTPUTS:
% --------
%   fxyn	  : (3*Np) 
%   zxy		  : (2,Np) [zx; zy]
%   DfxynDuxyn	  : (3,3,Np) [fx,x fx,y fx,z;
%  			    fy,x fy,y fy,z;
%  			    fz,x fz,y fz,z];

% error('This is work in progress and has to be debugged');

  Np = length(uxyn)/3;
  
  if nargin==6  % static or transient without initial state
      h  = 0;
      tp = 0;
      uxynp  = zeros(Np*3, 1);
      fxynp  = zeros(Np*3, 1);
      dfxynp = zeros(Np*3, Np*3);
      dfxynp(1:3:end, 1:3:end, :) = diag(Kt(1,:));
      dfxynp(1:3:end, 2:3:end, :) = diag(Kt(3,:));
      dfxynp(2:3:end, 1:3:end, :) = diag(Kt(3,:));
      dfxynp(2:3:end, 2:3:end, :) = diag(Kt(2,:));
      dfxynp(3:3:end, 3:3:end, :) = diag(kn);
  elseif nargin==11
      tp     = varargin{1};
      uxynp  = varargin{2};
      fxynp  = varargin{3};
      dfxynp = varargin{4};
      h      = varargin{5};
  else
      fprintf('%d inputs unknown\n',nargin);
      keyboard
      error(sprintf('%d inputs unknown',nargin));      
  end
  h = h(:);
  del_cst = [cos(h(h~=0)*t) sin(h(h~=0)*t)]'-[cos(h(h~=0)*tp) sin(h(h~=0)*tp)]';
  del_cst = [zeros(1, h(1)==0), del_cst(:)'];
  
  cst = [cos(h(h~=0)*t) sin(h(h~=0)*t)]';
  cst = [ones(1, h(1)==0), cst(:)'];
  
  if numel(Kt)==3 || length(mu)==1 || length(kn)==1 || length(gap)==1
      Kt  = repmat(Kt(:), 1, Np);
      kn  = kn*ones(1, Np);
      mu  = mu*ones(1, Np);
      gap = gap*ones(1, Np);
  end
  
  fxyn = zeros(Np*3, 1);
  DfxynDuxyn = zeros(Np*3, Np*3, length(del_cst));
  for di=1:Np
      xi = (di-1)*3+1;
      yi = (di-1)*3+2;
      ni = (di-1)*3+3;
      
      fxyn(ni) = max(kn(di)*(uxyn(ni)-gap(di)), 0);
      if fxyn(ni)==0  % separation
          fxyn(xi) = 0;
          fxyn(yi) = 0;
          
          DfxynDuxyn(ni, [xi yi ni], :) = 0;
          
          DfxynDuxyn(xi, [xi yi ni], :) = 0;
          
          DfxynDuxyn(yi, [xi yi ni], :) = 0;
      else  % contact
          
          % Update Normal Jacobian
          DfxynDuxyn(ni, [xi yi], :) = 0;
          DfxynDuxyn(ni, ni, :) = kn(di)*cst;
          
          % Predictor-Corrector for Tangential Forces
          Kmat = [Kt(1, di), Kt(3, di);
                  Kt(3, di), Kt(2, di)];
          fslip = mu(di)*fxyn(ni);
          dfSlipdn = mu(di)*DfxynDuxyn(ni, ni, :);

          % Stick Prediction
          fxystick = Kmat*(uxyn([xi;yi])-uxynp([xi;yi])) + fxynp([xi;yi]);
          fT = sqrt(fxystick'*fxystick);
          
          % Stick stiffness
          DfxyDuxynstick = [Kmat [0;0]].*permute(del_cst, [1, 3, 2]) + ...
              dfxynp([xi; yi], [xi yi ni], :);
          
          % Slip Correction
          if fT<fslip  % stick
              fxyn([xi;yi]) = fxystick;
              
              DfxynDuxyn([xi; yi], [xi yi ni], :) = DfxyDuxynstick;
          else
%               fxyn([xi;yi]) = fxystick;
              
%               DfxynDuxyn([xi; yi], [xi yi ni], :) = DfxyDuxynstick;
              
              fxyn([xi;yi]) = (fslip/fT)*fxystick;
              
              dfT = (fxystick(1)*DfxyDuxynstick(1,:,:) + ...
                     fxystick(2)*DfxyDuxynstick(2,:,:))/(fT);
              
              DfxynDuxyn([xi; yi], [xi yi ni], :) = (fslip/fT)*DfxyDuxynstick - ...
                  [(fslip/fT^2)*dfT*fxystick(1); ...
                   (fslip/fT^2)*dfT*fxystick(2)];

              DfxynDuxyn([xi; yi], ni, :) = DfxynDuxyn([xi; yi], ni, :) + ...
                  (dfSlipdn/fT).*fxystick;
          end
      end
  end
end
