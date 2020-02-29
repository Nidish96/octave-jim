function [T,X,Xd,Xdd] = HHTA_NONLIN(M,C,K,FN,Fnl,X0,Xd0,t0,t1,dt,a,b,g, ...
                                    opts)
%HHTA_NONLIN returns the HHT-Alpha time march solution for problems
%with nonlinearities
% USAGE:
%   [T,X,Xd] = HHTA(M,C,K,@(t) FN(t),X0,Xd0,t0,t1,h,a,g,b);
% INPUTS:
%   M       : NxN Intertia matrix
%   C       : NxN Proportional damping matrix
%   K       : NxN Stiffness matrix
%   FN      : Function handle for nonhomogeneous part: F(t),  returning
% 		Fex
%               Nx1 vector
%   Fnl     : Function handle for nonlinearity: F(t, X, Xd), returning 
% 		Fnl       , dFnldX    , dFnldXdot
%		Nx1 vector, NxN matrix, NxN matrix
%   X0      : Nx1 Vector of initial displacements
%   Xd0     : Nx1 Vector of initial velocities
%   t0      : Initial time
%   t1      : Final time
%   dt      : Time step
%   a       : Alpha parameter
%   b       : Beta parameter
%   g       : Gamma parameter
%   opts    : options structure with,
% 	reletol (float)
% 	etol    (float)
% 	utol    (float)
% 	rtol    (float)
% 	Display (boolean)
% OUTPUTS:
%   T       : 1xNt Time vector
%   X       : NxNt Displacement time series
%   Xd      : NxNt Velocity time series
%   Xdd     : NxNt Acceleration time series

    N = length(X0);
    X0 = reshape(X0,N,1);
    Xd0 = reshape(Xd0,N,1);
    
    Z1 = M+(1+a)*g*dt*C+(1+a)*b*dt^2*K;
    Z2 = M-(1+a)*(1-g)*dt*C-(1+a)*(0.5-b)*dt^2*K;
    Z3 = (1+a)*dt*K;
    
    T = t0:dt:t1;       Nt = length(T);
    X = zeros(N,Nt);    Xd = X;             Xdd = X;
    X(:,1) = X0;        Xd(:,1) = Xd0;      Xdd(:,1) = M\(FN(t0)-C*Xd0-K*X0-Fnl(t0,X0,Xd0));
    
    Xdd0 = Xdd(:, 1);
    Fnl_p = Fnl(T(1), X(:, 1), Xd(:, 1));  % Non-linear force at First step

    for i=2:length(T)
        % Explicit Predictor
        Xdd(:, i) = Xdd(:, i-1);
        
        % Corrector Iterations
        [Fnl_np1pa, dFdX, dFdXdot] = ...
            Fnl(T(i-1)+(1+a)*dt, ...
                X(:, i-1) + (1+a)*dt*Xd(:, i-1) + (1+a)*dt^2*((.5-b)*Xdd(:, i-1)+b*Xdd(:, i)), ...
                Xd(:, i-1) + (1+a)*dt^2*((1-g)*Xdd(:, i-1)+g*Xdd(:, i)));
        R = Z1*Xdd(:, i) - Z2*Xdd(:, i-1) + Z3*Xd(:, i-1) + ...
            (Fnl_np1pa-Fnl_p) - (FN(T(i)+(1+a)*dt)-FN(T(i-1)));  % Residue
        J = Z1 + (1+a)*(b*dt^2*dFdX + g*dt*dFdXdot);
        du = -J\R;
        e = abs(R'*du);
        r = mean(R.^2);
        u = mean(du.^2);
        
        e0 = e;
        r0 = r;
        u0 = u;
        it = 0;
        
        flag = 1*(e/e0<opts.reletol) + 2*(e<opts.etol) + 4*(r<opts.rtol) ...
               + 8*(u<opts.utol);
        if opts.Display
            fprintf('ITN, E, E/E0, r, du\n%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, flag);
        end        
        while (flag < 3) || (it==0) % 
            Xdd(:, i) = Xdd(:, i) + du;
            it = it+1;
            
            [Fnl_np1pa, dFdX, dFdXdot] = ...
                Fnl(T(i-1)+(1+a)*dt, ...
                    X(:, i-1) + (1+a)*dt*Xd(:, i-1) + (1+a)*dt^2*((.5-b)*Xdd(:, i-1)+b*Xdd(:, i)), ...
                    Xd(:, i-1) + (1+a)*dt^2*((1-g)*Xdd(:, i-1)+g*Xdd(:, i)));
            
            % Residual, Jacobian, and Update
            R = Z1*Xdd(:, i) - Z2*Xdd(:, i-1) + Z3*Xd(:, i-1) + ...
                (Fnl_np1pa-Fnl_p) - (FN(T(i-1)+(1+a)*dt)-FN(T(i-1)));  % Residue
            J = Z1 + (1+a)*(b*dt^2*dFdX + g*dt*dFdXdot);
            du = -J\R;
            e = abs(R'*du);
            r = mean(R.^2);
            u = mean(du.^2);
            
            flag = 1*(e/e0<opts.reletol) + 2*(e<opts.etol) + 4*(r<opts.rtol) ...
                   + 8*(u<opts.utol);
            if opts.Display
                fprintf('%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, flag);
            end
            
            if it>opts.ITMAX
                flag = 0;
                break;
            end
        end
        
        if flag == 0 || any(~isfinite(abs(X(:, i))))
            disp('No Convergence/Non finite march : Returning')

            X = X(:, 1:i-1);
            Xd = Xd(:, 1:i-1);
            Xdd = Xdd(:, 1:i-1);
            T = T(:, 1:i-1);
            
            break;
        end
        
        % Update Displacements and Velocities
        Xd(:,i) = Xd(:,i-1) + dt*((1-g)*Xdd(:,i-1)+g*Xdd(:,i));
        X(:,i) = X(:,i-1) + dt*Xd(:,i-1) + ...
                 dt^2*((0.5-b)*Xdd(:,i-1)+b*Xdd(:,i));
        Fnl_p = Fnl(T(i), X(:, i), Xd(:, i));  % Non-linear force
                                               % at step
        
        if opts.Display
            fprintf('---------------------------------------------------\n');
            fprintf('%.4e/%.4e %.4e\n', T(i), t1, dt);
            fprintf('---------------------------------------------------\n');
        end
    end
end