function [R, dRdU, dRdw] = UNIL_RESFUN(Uw, m, c, k, g, kg, Fl, Nt, h)
  Nhc = sum(h==0)+2*sum(h~=0);
  [E, dEdw] = HARMONICSTIFFNESS(m, c, k, Uw(end), h);
  cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);
  
  u = TIMESERIES_DERIV(Nt, h, Uw(1:end-1), 0);
  fnl = max(kg*(u-g), 0);  cont = (fnl~=0);
  jnl = kg.*cont;
  
  R = E*Uw(1:end-1) + GETFOURIERCOEFF(h, fnl) - Fl;
  dRdU = E + GETFOURIERCOEFF(h, jnl.*cst);
  dRdw = dEdw*Uw(1:end-1);
end
