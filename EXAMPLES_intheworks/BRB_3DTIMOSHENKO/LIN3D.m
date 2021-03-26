function [fnl, dfnldu] = LIN3D(t, u, ktx, kty, kn, mu, gap, varargin)
%ELDRYFRICT3D returns the force and jacobian for the 3D elastic dry
%friction element
%  
%  USAGE:
%    [fnl, dfnldu] = ELDRYFRICT3D(t, u, ktx, kty, kn, mu, varargin)
%  INPUT:
%       t       : scalar
%       u       : 3*Np x 1  [u1; v1; w1; u2; v2; w2; ...]
%       ktx     : scalar or Np x 1
%       kty     : scalar or Np x 1
%       kn      : scalar or Np x 1
%       mu      : scalar or Np x 1
%       gap	: scalar or Np x 1
%       h       : Nh x 1
%       tp      : scalar
%       up      : 3*Np x 1  [same as u]
%       fp      : 3*Np x 1  [fx1; fy1; fz1; fx2; fy2; fz2; ...]
%       dfp     : 3*Np x 3*Np x Nhc
%  OUTPUTs:
%       fnl     : 3*Np x 1 
%       dfnldu  : 3*Np x 3*Np x Nhc

  Np = length(u)/3;
  if Np~=fix(Np)
    fprintf('Non-integer number of points\n')
    keyboard
    error('Non-integer number of points')
  end
  
  if nargin==7  % possibly quasi-static operation
    h = 0; 
    tp = 0;
    up = zeros(1, 3*Np);
    fp = zeros(1, 3*Np);
    dfp = zeros(3*Np, 3*Np);
    dfp(sub2ind(size(dfp), 1:3:(3*Np), 1:3:(3*Np))) = ktx;
    dfp(sub2ind(size(dfp), 2:3:(3*Np), 2:3:(3*Np))) = kty;
    dfp(sub2ind(size(dfp), 3:3:(3*Np), 3:3:(3*Np))) = kn;
  elseif nargin==12
    h = varargin{1};
    tp = varargin{2};
    up = varargin{3};
    fp = varargin{4};
    dfp = varargin{5};
  else
    fprintf('%d inputs unknown\n', nargin);
    keyboard
    error(sprintf('%d inputs unknown', nargin));
  end
  u = u(:);
  up = up(:);
  fp = fp(:);
  
  h = h(:);  Nhc = sum(h==0)+2*sum(h~=0);
  del_cst = [cos(h(h~=0)*t) sin(h(h~=0)*t)]'-[cos(h(h~=0)*tp) sin(h(h~=0)*tp)]';
  del_cst = [zeros(1, h(1)==0), del_cst(:)'];  % 1xNhc
  del_cst = permute(del_cst, [1 3 2]);

  cst = [cos(h(h~=0)*t) sin(h(h~=0)*t)]';
  cst = [ones(1, h(1)==0), cst(:)'];  % 1xNhc
  cst = permute(cst, [1 3 2]);  
  
  if length(ktx)==1 || length(kty)==1 || length(kn)==1 || length(mu)==1 || length(gap)==1
    ktx = ktx(1)*ones(Np,1);
    kty = kty(1)*ones(Np,1);
    kn  = kn(1)*ones(Np,1);
    mu  = mu(1)*ones(Np,1);
    gap  = gap(1)*ones(Np,1);
  end

  fnl    = zeros(Np*3, 1);
  dfnldu = zeros(Np*3, Np*3, Nhc);
  
  fnl(1:3:end) = ktx.*(u(1:3:end)-up(1:3:end))+fp(1:3:end);
  fnl(2:3:end) = kty.*(u(2:3:end)-up(2:3:end))+fp(2:3:end);
  fnl(3:3:end) = kn.*u(3:3:end);
  
  dfnldu(1:3:end, 1:3:end, :) = diag(ktx).*del_cst+dfp(1:3:end, 1:3:end, :);
  dfnldu(2:3:end, 2:3:end, :) = diag(kty).*del_cst+dfp(2:3:end, 2:3:end, :);
  dfnldu(3:3:end, 3:3:end, :) = diag(kn).*cst;
  
  fnl(1:3:end) = ktx.*u(1:3:end);
  fnl(2:3:end) = kty.*u(2:3:end);
  fnl(3:3:end) = kn.*u(3:3:end);
  
  dfnldu(1:3:end, 1:3:end, :) = diag(ktx).*cst;
  dfnldu(2:3:end, 2:3:end, :) = diag(kty).*cst;
  dfnldu(3:3:end, 3:3:end, :) = diag(kn).*cst;  
  
%   if dfnldu(1:3:end, 1:3:end, 1) == 0
%       dfnldu(1:3:end, 1:3:end, 1) = diag(ktx);
%   end
%   if dfnldu(2:3:end, 2:3:end, 1) == 0
%       dfnldu(2:3:end, 2:3:end, 1) = diag(kty);
%   end
end