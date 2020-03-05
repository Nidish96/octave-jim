function [fxyn, zxy, DfxynDuxyn, DfxynDuxynd, DfxynDkxynmu] = ELDRYFRICT(uxyn, zxy, uxynd, kxynmu, N0)
%ELDRYFICT returns the forces and jacobians for the elastic dry friction model
%
% USAGE:
% ------  
%   [fxyn, zxy, DfxynDuxyn, DfxynDuxynd] = ELDRYFRICT(uxyn, zxy, uxynd, kxynmu, N0);
% INPUTS:
% -------
%   uxyn	  : (3,Np) [ux; uy; un]
%   uxynd	  : (3,Np) [uxd; uyd; und]  
%   zxy		  : (2,Np) [zx; zy]
%   kxynmu	  : (4,Np) [kx; ky; kn; mu]
%   N0		  : (1,Np) [N0]
% OUTPUTS:
% --------
%   fxyn	  : (3,Np) [fx; fy; fn]
%   zxy		  : (2,Np) [zx; zy]
%   DfxynDuxyn	  : (3,3,Np) [fx,x fx,y fx,z;
%  			    fy,x fy,y fy,z;
%  			    fz,x fz,y fz,z];
%   DfxynDuxynd   : (3,3,Np) [fx,xd fx,yd fx,zd;
%  			     fy,xd fy,yd fy,zd;
%  			     fz,xd fz,yd fz,zd];
%   DfxynDkxynmu  : (3,4,Np) [fx,kx fx,ky fx,kn fx,mu;
%                           fy,kx fy,ky fy,kn fy,mu;
%                           fn,kx fn,ky fn,kn fn,mu];

  Np = size(uxyn, 2);
  fxyn = zeros(3, Np);
  DfxynDuxyn = zeros(3, 3, Np);
  DfxynDuxynd = zeros(3, 3, Np);
  DfxynDkxynmu = zeros(3, 4, Np);
  
				% 1. STICK (PREDICTION)
  fxyn(1, :) = kxynmu(1,:).*(uxyn(1, :)-zxy(1, :));
  fxyn(2, :) = kxynmu(2,:).*(uxyn(2, :)-zxy(2, :));
  fxyn(3, :) = kxynmu(3,:).*max(uxyn(3, :)+N0./kxynmu(3,:), 0);
					      % Derivatives
  DfxynDuxyn(1, 1, :) = kxynmu(1,:); % fx,ux
  DfxynDuxyn(2, 2, :) = kxynmu(2,:); % fx,uy
  DfxynDuxyn(3, 3, :) = kxynmu(3,:); % fx,un
				     % Parameter derivatives
  DfxynDkxynmu(1, 1, :) = uxyn(1,:)-zxy(1,:); % fx,kx
  DfxynDkxynmu(2, 2, :) = uxyn(2,:)-zxy(2,:); % fy,ky
  DfxynDkxynmu(3, 3, :) = uxyn(3,:); % fn,kn

				% 2. SEPARATION
  isep = find(fxyn(3,:)==0);  % indices of separated points
  zxy(:, isep) = uxyn(1:2, isep);
				% Everything is zero when separated
  fxyn(:, isep) = 0;
  DfxynDuxyn(:, :, isep) = 0;
  DfxynDuxynd(:, :, isep) = 0;
  DfxynDkxynmu(:, :, isep) = 0;

				      % 3. SLIP
  fT = sqrt(sum(fxyn(1:2, :).^2, 1));  % Tangential force magnitude
  fslip = kxynmu(4,:).*fxyn(3, :);        % Slip force magnitude

  islips = find(fT>fslip);       % indices of slipped points
  islips = setdiff(islips, isep);% Not interested in separated points - everything's zero there
  
  fT(fT<eps) = 1.0;  % Avoid dividing by zeros
  if length(islips) ~= 0
				% Derivatives
    DfxynDuxyn(1, 1, islips) = fslip(islips).*fxyn(2,islips).^2.*kxynmu(1,islips)./fT(islips).^3;% fx,ux
    DfxynDuxyn(1, 2, islips) = -fslip(islips).*prod(fxyn(1:2,islips),1).*kxynmu(2,islips)./fT(islips).^3;% fx,uy
    DfxynDuxyn(1, 3, islips) = kxynmu(4,islips).*kxynmu(3,islips).*fxyn(1,islips)./fT(islips);% fx,un
    
    DfxynDuxyn(2, 1, islips) = -fslip(islips).*prod(fxyn(1:2,islips),1).*kxynmu(1,islips)./fT(islips).^3;% fy,ux
    DfxynDuxyn(2, 2, islips) = fslip(islips).*fxyn(1,islips).^2.*kxynmu(2,islips)./fT(islips).^3;% fy,uy
    DfxynDuxyn(2, 3, islips) = kxynmu(4,islips).*kxynmu(3,islips).*fxyn(2,islips)./fT(islips);% fx,un
				% Parameter Derivatives
    DfxynDkxynmu(1, 1, islips) = fslip(islips).*fxyn(2,islips).^2.*(uxyn(1,islips)-zxy(1,islips))./fT(islips).^3;% fx,kx
    DfxynDkxynmu(1, 2, islips) = -fslip(islips).*prod(fxyn(1:2,islips),1).*(uxyn(2,islips)-zxy(2,islips))./fT(islips).^3;% fx,ky
    DfxynDkxynmu(1, 3, islips) = kxynmu(4,islips).*uxyn(3,islips).*fxyn(1,islips)./fT(islips);% fx,kn
    DfxynDkxynmu(1, 4, islips) = (fslip(islips)./kxynmu(4,islips)).*fxyn(1,islips)./fT(islips);% fx,mu

    DfxynDkxynmu(2, 1, islips) = -fslip(islips).*prod(fxyn(1:2,islips),1).*(uxyn(1,islips)-zxy(1,islips))./fT(islips).^3;% fy,kx  
    DfxynDkxynmu(2, 2, islips) = fslip(islips).*fxyn(1,islips).^2.*(uxyn(2,islips)-zxy(2,islips))./fT(islips).^3;% fy,ky
    DfxynDkxynmu(2, 3, islips) = kxynmu(4,islips).*uxyn(3,islips).*fxyn(2,islips)./fT(islips);% fy,kn
    DfxynDkxynmu(2, 4, islips) = (fslip(islips)./kxynmu(4,islips)).*fxyn(2,islips)./fT(islips);% fy,mu

				% Update forces to lie on slip-cone
    fxyn(1:2, islips) = fxyn(1:2, islips).*repmat(fslip(islips)./fT(islips), 2, 1);

    zxy(:, islips) = uxyn(1:2, islips) - fxyn(1:2, islips)./kxynmu(1:2, islips);
  end

  % Remove N0 (unbalanced static force, assumed to just exist to activate the sliders)
  fxyn(3, :) = fxyn(3, :) - N0;
  fxyn(3, isep) = 0;
end
