function [fnl, dfnldu] = ELDRYFRICT2D(t, utun, kt, kn, mu, gap, varargin)
%ELDRYFRICT2D returns the force and jacobian for the Jenkins element
%
%  USAGE:
%    [fnl, dfnldu] = ELDRYFRICT2D(t, utun, kt, kn, mu, gap, h, tp, up, fp, dfp)
%  INPUT:
%       t       : scalar
%       utun    : 2*Nd x 1
%       kt      : scalar or Nd x 1
%       kn      : scalar or Nd x 1
%       mu      : scalar or Nd x 1
%       gap     : scalar or Nd x 1
%       h       : Nh x 1
%       tp      : scalar
%       up      : 2*Nd x 1
%       fp      : 2*Nd x 1
%       dfp     : 2*Nd x 2*Nd x Nhc
%  OUTPUTs:
%       fnl     : 2*Nd x 1 
%       dfnldu  : 2*Nd x 2*Nd x Nhc

    Nd = length(utun)/2;
    if nargin==6  % static or transient operatin without initial state
      h = 0;
      tp = 0;
      up = zeros(Nd*2, 1);
      fp = zeros(Nd*2, 1);
      dfp = permute(kron(eye(Nd), diag([kt, kn])), [3 1 2]);
    elseif nargin==11
      h = varargin{1};
      tp = varargin{2};
      up = varargin{3};
      fp = varargin{4};
      dfp = varargin{5};
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
    
    if length(kt)==1 || length(mu)==1 || length(kn)==1 || length(gap)==1
      kt  = kt*ones(Nd, 1);
      kn  = kn*ones(Nd, 1);
      mu  = mu*ones(Nd, 1);
      gap = gap*ones(Nd, 1);
    end

    fnl = zeros(Nd*2, 1);
    dfnldu = zeros(Nd*2, Nd*2, length(del_cst));
    for di=1:Nd  % Iterate over the 
        ti = (di-1)*2 + 1;  % Tangential DOf ID
        ni = (di-1)*2 + 2;  % Normal DOF ID
        
        fnl(ni) = max(kn(di)*(utun(ni)-gap(di)), 0);
        if fnl(ni)==0  % Separation
            fnl(ti) = 0;
            dfnldu(ni, ni, :) = 0;
            dfnldu(ti, ti, :) = 0;
            dfnldu(ti, ni, :) = 0;
        else  % Contact
            dfnldu(ni, ni, :) = kn(di).*cst;  % Direct formula (check which is better)
            
            fslip = mu(di)*fnl(ni);  % >=0 because of unilateral contact
            fsp = kt(di)*(utun(ti)-up(ti))+fp(ti);  % Stick prediction
            
            if abs(fsp)<fslip  % Stick
                fnl(ti) = fsp;
                dfnldu(ti, ti, :) = kt(di).*del_cst+squeeze(dfp(1, ti, ti, :))';
                dfnldu(ti, ni, :) = squeeze(dfp(1, ti, ni, :));
            else  % Slip
                fnl(ti) = fslip*sign(fsp);
                
                dfnldu(ti, ti, :) = 0;
                dfnldu(ti, ni, :) = (mu(di)*kn(di)*sign(fsp)).*cst;
            end            
        end
    end
end