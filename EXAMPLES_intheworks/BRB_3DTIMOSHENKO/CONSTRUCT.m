function [M, K] = CONSTRUCT(Beam, pars)
  Ne = length(Beam.X)-1;

  M = zeros((Ne+1)*6);
  K = zeros((Ne+1)*6);
  for e = 1:Ne
    if Beam.X(e+1)-Beam.X(e)<0
      error('Negative Jacobian - Quitting')
    end
    [Me, Ke] = TM3D_ELMATS(Beam.X(e+1)-Beam.X(e), pars);

    is = (e-1)*6+1;
    ie = (e+1)*6;
    M(is:ie, is:ie) = M(is:ie, is:ie) + Me;
    K(is:ie, is:ie) = K(is:ie, is:ie) + Ke;
  end
end
