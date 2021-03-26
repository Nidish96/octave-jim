function [fnl, dfnldu] = ELDRYFRICT3D(t, u, ktx, kty, kn, mu, gap, varargin)
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
%       tp      : scalar (previous time instant)
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
    up = zeros(3*Np);
    fp = zeros(3*Np);
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
  for pi=1:Np
    ux  = u((pi-1)*3+1);   uy  = u((pi-1)*3+2);   uz  = u((pi-1)*3+3);
    uxp = up((pi-1)*3+1);  uyp = up((pi-1)*3+2);  uzp = up((pi-1)*3+3);

    fxp = fp((pi-1)*3+1);  fyp = fp((pi-1)*3+2);  fzp = fp((pi-1)*3+3);

    dfxxp = dfp((pi-1)*3+1, (pi-1)*3+1, :);
    dfxyp = dfp((pi-1)*3+1, (pi-1)*3+2, :);
    dfxzp = dfp((pi-1)*3+1, (pi-1)*3+3, :);
    
    dfyxp = dfp((pi-1)*3+2, (pi-1)*3+1, :);
    dfyyp = dfp((pi-1)*3+2, (pi-1)*3+2, :);
    dfyzp = dfp((pi-1)*3+2, (pi-1)*3+3, :);
    
    dfzzp = dfp((pi-1)*3+3, (pi-1)*3+3, :);

    fz = max(kn(pi)*(uz-gap(pi)), 0);
    dfzz = (fz>0)*kn(pi)*cst;
%     
%     fz = 1e6;
%     dfzz = dfzz*0;

    if fz==0  % separation
      fx = 0;  fy = 0;

      dfxx = 0; dfxy = 0; dfxz = 0;
      dfyx = 0; dfyy = 0; dfyz = 0;
    else  % contact
      % stick-prediction
      fx_pred = ktx(pi)*(ux-uxp)+fxp;
      fy_pred = kty(pi)*(uy-uyp)+fyp;
      fT_pred = sqrt(fx_pred^2+fy_pred^2);

      dfxx_pred = ktx(pi)*del_cst + dfxxp;
      dfxy_pred = dfxyp;
      dfxz_pred = dfxzp;

      dfyx_pred = dfyxp;
      dfyy_pred = kty(pi)*del_cst + dfyyp;
      dfyz_pred = dfyzp;
      
      % slip limit
      fslip = mu(pi)*fz;

      if fT_pred<fslip  % stick
        fx = fx_pred;  fy = fy_pred;

        dfxx = dfxx_pred;
        dfxy = dfxy_pred;
        dfxz = dfxz_pred;

        dfyx = dfyx_pred;
        dfyy = dfyy_pred;
        dfyz = dfyz_pred;
      else  % slip
%           fprintf('%d slips!\n', pi)
        fx = fslip*fx_pred/fT_pred;
        fy = fslip*fy_pred/fT_pred;

        dfxx = fslip*fy_pred*(fy_pred*dfxx_pred - fx_pred*dfyx_pred)/fT_pred^3;
        dfxy = fslip*fy_pred*(fy_pred*dfxy_pred - fx_pred*dfyy_pred)/fT_pred^3;
%         dfxz = mu(pi)*fx_pred/fT_pred*dfzz;
        dfxz = mu(pi)*fx_pred/fT_pred*dfzz + fslip*fy_pred*(fy_pred*dfxz_pred - fx_pred*dfyz_pred)/fT_pred^3;

        dfyx = fslip*fx_pred*(fx_pred*dfyx_pred - fy_pred*dfxx_pred)/fT_pred^3;
        dfyy = fslip*fx_pred*(fx_pred*dfyy_pred - fy_pred*dfxy_pred)/fT_pred^3;
%         dfyz = mu(pi)*fy_pred/fT_pred*dfzz;
        dfyz = mu(pi)*fy_pred/fT_pred*dfzz + fslip*fx_pred*(fx_pred*dfyz_pred - fy_pred*dfxz_pred)/fT_pred^3;
      end
    end
    
    % assembly
    if ~isfinite(fx) || ~isfinite(fy) || ~isfinite(fz)
        keyboard
    end
    fnl((pi-1)*3+(1:3)) = [fx; fy; fz];

    dfnldu((pi-1)*3+1, (pi-1)*3+1, :) = dfxx;
    dfnldu((pi-1)*3+1, (pi-1)*3+2, :) = dfxy;
    dfnldu((pi-1)*3+1, (pi-1)*3+3, :) = dfxz;

    dfnldu((pi-1)*3+2, (pi-1)*3+1, :) = dfyx;
    dfnldu((pi-1)*3+2, (pi-1)*3+2, :) = dfyy;
    dfnldu((pi-1)*3+2, (pi-1)*3+3, :) = dfyz;

    dfnldu((pi-1)*3+3, (pi-1)*3+3, :) = dfzz;
  end
end

