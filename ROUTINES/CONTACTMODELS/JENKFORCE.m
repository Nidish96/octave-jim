function [fnl, dfnldu] = JENKFORCE(t, u, kt, muN, varargin)
%JENKFORCE returns the force and jacobian for the Jenkins element
%
%  USAGE:
%    [fnl, dfnldu] = JENKFORCE(t, u, kt, muN, h, tp, up, fp, dfp)
%  INPUT:
%       t       : scalar
%       u       : Nd x 1
%       kt      : scalar or Nd x 1
%       muN     : scalar or Nd x 1
%       h       : Nh x 1
%       tp      : scalar
%       up      : Nd x 1
%       fp      : Nd x 1
%       dfp     : Nd x Nd x Nhc
%  OUTPUTs:
%       fnl     : Nd x 1 
%       dfnldu  : Nd x Nd x Nhc
    
    if nargin==4
      h = 0;
      tp = 0;
      up = zeros(size(u));
      fp = zeros(size(u));
      dfp = permute(kt.*eye(length(u)), [3 1 2]);
    elseif nargin==9
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
    
    if length(kt)==1 && length(muN) == 1
      kt = kt*ones(length(u), 1);
      muN = muN*ones(length(u), 1);
    end

    fnl    = zeros(length(u));
    dfnldu = zeros(length(u), length(u), length(del_cst));
    for di=1:length(u)
      fsp = kt(di)*(u(di)-up(di))+fp(di);   % Stick prediction
      
      fnl(di) = fsp*(abs(fsp)<muN(di)) + muN(di)*sign(fsp)*(abs(fsp)>=muN);  % Nonlinear force
      
      dfnldu(di, di, :) = (kt(di).*del_cst+squeeze(dfp(1, di, di, :))').*(abs(fsp)<muN(di));  % Nonlinear Jacobian
    end
end