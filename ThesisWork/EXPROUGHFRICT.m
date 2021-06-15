function [fxyn, varargout] = EXPROUGHFRICT(t, uxyn, ct, cn, lam, mu, gap, varargin)
%EXPROUGHFRICT returns the forces and jacobians for the elastic dry friction model
%
% USAGE:
% ------  
%   [fxyn, zxy, DfxynDuxyn] = ...
%      EXPROUGHFRICT(t, uxyn, Kt, kn, mu, gap, tp, uxynp, fxynp, dfxynp, h)
% INPUTS:
% -------
%   uxyn	  : (3*Np, 1) [[ux; uy; un]_1; [ux; uy; un]_2; ...]
%   uxynd	  : (3*Np, 1) [[uxd; uyd; und]_1; [uxd; uyd; und]_2; ...]
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

%     fxyn = 0;
%     DfxynDuxyn = 0;
%     return;

  Np = length(uxyn)/3;
  
  if nargin==7  % static or transient without initial state
      h  = 0;
      tp = 0;
      uxynp  = zeros(Np*3, 1);
      fxynp  = zeros(Np*3, 1);
      dfxynp = zeros(Np*3, Np*3);
  elseif nargin==12
      tp     = varargin{1};
      uxynp  = varargin{2}(:);
      fxynp  = varargin{3}(:);
      dfxynp = squeeze(varargin{4});
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
  
  if numel(ct)==3 || numel(cn)==3 || length(mu)==1 || length(lam)==1 || length(gap)==1
      cn  = cn*ones(1, Np);
      ct  = ct*ones(1, Np);
      lam = lam*ones(1, Np);
      mu  = mu*ones(1, Np);
      gap = gap*ones(1, Np);
  end
  
  fxyn = zeros(Np*3, 1);
  DfxynDuxyn = zeros(Np*3, Np*3, length(del_cst));
  
  fxyn(3:3:end) = cn.*exp(lam.*(uxyn(3:3:end)-gap));
  if length(del_cst)==1
    DfxynDuxyn(sub2ind([Np*3 Np*3], 3:3:Np*3, 3:3:Np*3)) = lam.*fxyn(3:3:end);
  end
  if nargout==1  % Only force required (can be vectorized easily)
      % Stick Prediction
      kt = ct.*exp(lam.*(uxyn(3:3:end)-gap));
      fxyn(1:3:end) = kt.*(uxyn(1:3:end)-uxynp(1:3:end))+fxynp(1:3:end);
      fxyn(2:3:end) = kt.*(uxyn(2:3:end)-uxynp(2:3:end))+fxynp(2:3:end);
      
      fT = sqrt(fxyn(1:3:end).^2+fxyn(2:3:end).^2);
      fslips = mu.*fxyn(3:3:end);
      
      sis = find(fT>fslips);
      xsis = (sis-1)*3+1;  ysis = (sis-1)*3+2;
      
      fxyn(xsis) = fslips(sis).*fxyn(xsis)./fT(sis);
      fxyn(ysis) = fslips(sis).*fxyn(ysis)./fT(sis);
  else  % Force with Jacobian required
      for di=1:Np  % Can Parallelize or at least vectorize for speed
          xi = (di-1)*3+1;
          yi = (di-1)*3+2;
          ni = (di-1)*3+3;

          % Normal Contact
          if length(del_cst)>1
            DfxynDuxyn(ni, ni, :) = lam(di)*fxyn(ni)*cst;
          end

          % Tangential Stiffness
          kt = ct(di)*exp(lam(di)*(uxyn(ni)-gap(di)));
          dktdun = lam(di)*kt;

          % Predictor-Corrector for Tangential Forces
          fslip = mu(di)*fxyn(ni);
          dfSlipdn = lam(di)*fslip*cst;  %(4)

          % Stick Prediction
          del_u = uxyn([xi;yi])-uxynp([xi;yi]);
          fxystick = kt*del_u + fxynp([xi;yi]);  %(3)
          fT = sqrt(fxystick'*fxystick);

          % Stick stiffness
          % zeroth harmonic gradient set to instantaneous gradient $(1)
          if length(del_cst)>1
              DfxyDuxynstick = [eye(2)*kt [0;0]].*permute(del_cst, [1, 3, 2]) + ...
                  [zeros(2) dktdun*del_u].*permute(cst, [1, 3, 2]) + ...
                  dfxynp(xi:yi, xi:ni, :);
              DfxyDuxynstick(:, :, 1) = [eye(2)*kt dktdun*del_u]; 
          else
              DfxyDuxynstick = [[kt 0;0 kt] dktdun*del_u]; 
          end

          % Slip Correction
          if fT<fslip || fT<eps % stick
              fxyn([xi;yi]) = fxystick;

              DfxynDuxyn(xi:yi, xi:ni, :) = DfxyDuxynstick;  %(2)
          else
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
      varargout{1} = DfxynDuxyn;
  end
  
end