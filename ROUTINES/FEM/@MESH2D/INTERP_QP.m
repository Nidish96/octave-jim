function Uqp = INTERP_QP(m, U)
%INTERP_QP interpolates the given nodal data to quadrature points
  Uqp = m.Qm*U;
end
