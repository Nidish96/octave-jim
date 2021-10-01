function [fnl, dfnldu] = ELDRYFRICT2D(t, u, kt, kn, mu, gap, varargin)
%ELDRYFRICT2D returns the force and jacobian for the Jenkins element
%
%  USAGE:
%    [fnl, dfnldu] = ELDRYFRICT2D(t, u, kt, kn, mu, gap, h, tp, up, fp, dfp)
%  INPUT:
%       t       : scalar
%       u      : Nd x 1
%       kt      : scalar or Nd x 1
%       kn      : scalar or Nd x 1
%       mu      : scalar or Nd x 1
%       h       : Nh x 1
%       tp      : scalar
%       up      : Nd x 1
%       fp      : Nd x 1
%       dfp     : Nd x Nd x Nhc
%  OUTPUTs:
%       fnl     : Nd x 1 
%       dfnldu  : Nd x Nd x Nhc
    
    if nargin==6
      h = 0;
      tp = 0;
      up = zeros(size(u));
      fp = zeros(size(u));
      dfp = permute(kt.*eye(length(u)), [3 1 2]);
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
    
    u = u(:);
    ux = u(1:2:end);
    un = u(2:2:end);
    Nd = length(ux);
    
    h = h(:);
    hnnz = h(h~=0);
    
    del_cst = [cos(hnnz*t) sin(hnnz*t)]'-[cos(hnnz*tp) sin(hnnz*tp)]';
    del_cst = [zeros(1, h(1)==0), del_cst(:)'];
    
    cst = [cos(hnnz*t) sin(hnnz*t)]';
    cst = [ones(1, h(1)==0), cst(:)'];  % 1xNhc
    cst = permute(cst, [1 3 2]);  
    
    if length(kt)==1 && length(mu) == 1 && length(kn) == 1 && length(gap) == 1
      kt  = kt*ones(Nd, 1);
      kn  = kn*ones(Nd, 1);
      mu  = mu*ones(Nd, 1);
      gap = gap*ones(Nd, 1);
    end

    dfp = squeeze(dfp);
    
    fnl = zeros(2*Nd, 1);
    dfnldu = zeros(2*Nd, 2*Nd, size(dfp, 3));
    
    % Normal Contact
    fnl(2:2:end) = max(kn.*(un-gap), 0);
    dfnldu(2:2:end, 2:2:end, :) = diag(kn.*(fnl(2:2:end)~=0)).*cst;
    
    
    nt = zeros(Nd,1);  % 0-separation 1-contact 
    nt = (fnl(2:2:end)~=0);
    
    st = zeros(Nd,1);  % 0-stick 1-slip 
    for di=1:Nd
      xi = (di-1)*2+1;
      ni = (di-1)*2+2;
      
      if fnl(ni)~=0  % in contact
        fsp = kt(di)*(u(xi)-up(xi))+fp(xi);   % Stick prediction
        fslip = abs(mu(di)*fnl(ni));
        
        if abs(fsp)<fslip  % sticking
            st(di) = 0;
            
            fnl(xi) = fsp;
            
            dfnldu(xi, xi, :) = kt(di).*del_cst+squeeze(dfp(xi, xi, :))';
            dfnldu(xi, ni, :) = squeeze(dfp(xi, ni, :))';
        elseif abs(fsp)>=fslip  % slipping
            st(di) = 1;
            
            fnl(xi) = fslip*sign(fsp);
            
            dfnldu(xi, xi, :) = 0;
            dfnldu(xi, ni, :) = mu(di)*kn(di)*sign(fsp).*cst;
        end
      
%         st(di) = (abs(fsp)>fslip);
%         fnl(xi) = fsp*(abs(fsp)<=fslip) + fslip*sign(fsp)*(abs(fsp)>fslip);  % Nonlinear force
%       
%         dfnldu(xi, xi, :) = (kt(di).*del_cst+squeeze(dfp(xi, xi, :))').*(abs(fsp)<=fslip);  % Nonlinear Jacobian
%         dfnldu(xi, ni, :) = (squeeze(dfp(xi, ni, :))').*(abs(fsp)<=fslip) + ...
%                             (kn(di).*cst*sign(fsp)).*(abs(fsp)>fslip);  % Nonlinear Jacobian
      end
    end
%     if ~isempty(find(nt==0))
%       keyboard 
%     end
%     
%     disp([nt'; st'])
%     fprintf('\n');
%     
%     fnl(2:2:end) = kn.*un;
%     dfnldu(2:2:end, 2:2:end, :) = diag(kn).*cst;
% 
%     fnl(1:2:end) = kt.*ux; 
%     dfnldu(1:2:end, 1:2:end, :) = diag(kt).*cst;
end