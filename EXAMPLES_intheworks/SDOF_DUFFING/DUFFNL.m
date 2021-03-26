function [f, dfdx, dfdxd] = DUFFNL(t, x, xd, al)
  f = al*x.^3;
  dfdx = 3*al*x.^2;
  dfdxd = x*0;
end