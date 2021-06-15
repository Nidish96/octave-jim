function [fxyn, DfxynDuxyn] = EXPROUGHFRICT_sca(t, ux, uy, un, ct, cn, lam, mu, gap)
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
%   if nargin==9  % static or transient without initial state
      h  = 0;
      tp = 0;
      uxp = 0; uyp = 0; unp = 0;
      fxp = 0; fyp = 0; fnp = 0;
      
      dfxynp = [0 0 0;0 0 0;0 0 0];
%   elseif nargin==14
%       tp     = varargin{1};
%       uxynp  = [varargin{2}; varargin{3}; varargin{4}];
%       fxynp  = [varargin{5}; varargin{6}; varargin{7}];
%       dfxynp = squeeze(varargin{8});
%       h      = varargin{9};
%   else
%       fprintf('%d inputs unknown\n',nargin);
%       keyboard
%       error(sprintf('%d inputs unknown',nargin));      
%   end

  h = h(:);
  del_cst = [cos(h(h~=0)*t) sin(h(h~=0)*t)]'-[cos(h(h~=0)*tp) sin(h(h~=0)*tp)]';
  del_cst = [zeros(1, h(1)==0), del_cst(:)'];
  
  cst = [cos(h(h~=0)*t) sin(h(h~=0)*t)]';
  cst = [ones(1, h(1)==0), cst(:)'];
  
  fxyn = zeros(3, 1);
  DfxynDuxyn = zeros(3, 3, length(del_cst));
  uxyn = [ux; uy; un];
  
      xi = 1;
      yi = 2;
      ni = 3;
      
      % Normal Contact
      fxyn(ni) = cn*exp(lam*(uxyn(ni)-gap));
      DfxynDuxyn(ni, ni, :) = lam*fxyn(ni)*cst;
      
      % Tangential Stiffness
      kt = ct*exp(lam*(uxyn(ni)-gap));
      dktdun = lam*kt;
          
      % Predictor-Corrector for Tangential Forces
      Kmat = eye(2)*kt;
      dKmatdun = eye(2)*dktdun;
      fslip = mu*fxyn(ni);
      dfSlipdn = mu*DfxynDuxyn(ni, ni, :);

      % Stick Prediction
      fxystick = Kmat*(uxyn([xi;yi])-[uxp;uyp]) + [fxp;fyp];
      fT = sqrt(fxystick'*fxystick);
          
      % Stick stiffness
      DfxyDuxynstick = [Kmat [0;0]].*permute(del_cst, [1, 3, 2]) + ...
          [zeros(2) dKmatdun*(uxyn([xi;yi])-[uxp;uyp])].*permute(cst, [1, 3, 2]) + ...
          dfxynp([xi; yi], [xi yi ni], :);

      % zeroth harmonic gradient set to instantaneous gradient
      DfxyDuxynstick(:, :, 1) = [Kmat dKmatdun*(uxyn([xi;yi])-[uxp;uyp])]; 
          
      % Slip Correction
      if fT<fslip || fT<eps % stick
          fxyn([xi;yi]) = fxystick;

          DfxynDuxyn([xi; yi], [xi yi ni], :) = DfxyDuxynstick;
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
      
      fxyn = fxyn(1);
end
