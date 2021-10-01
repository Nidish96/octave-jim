function [R, dRdq0w2, dRdQ] = RQHBRESFUN(q0w2Q, fnl, jnl, Nt)
%RQHBRESFUN returns residual for RQNMA HB averaging (no rate-dependent
%nonlinearities
%
%  USAGE:
%   [R, dRdq0wn, dRdQ] = RQHBRESFUN(q0wnQ, fnl, Nt, h)
%  INPUTS:
%   q0wnQ   : 3x1 vector of [q0; wn; Q]
%   fnl     : fnl(t, q) force function handle returning Ntx1 vector
%   jnl     : jnl(t, q) jacobian function handle returning Ntx1 vector
%   Nt      : Time samples for AFT

    tau = linspace(0, 2*pi, Nt+1)'; tau(end) = [];
    
    q0 = q0w2Q(1);
    w2 = q0w2Q(2);
    w = sqrt(w2);
    Q = q0w2Q(3);
    
    qt = q0 + Q*cos(tau);
    cst = [ones(Nt,1) cos(tau)];
    
    ft = fnl(w*tau, qt);
    jt = jnl(w*tau, qt);
    
    Fh = GETFOURIERCOEFF(0:1, ft);
    Jh = GETFOURIERCOEFF(0:1, jt.*cst); 
    
    R = [0; -w2*Q] + Fh(1:2);
    dRdq0w2 = [0 0; 0 -Q] + [Jh(1:2,1), zeros(2,1)];
    dRdQ = [0; -w2] + Jh(1:2, 2);
end