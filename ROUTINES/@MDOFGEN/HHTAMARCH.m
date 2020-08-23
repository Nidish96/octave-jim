function [T, U, Ud, Udd, m] = HHTAMARCH(m, T0, T1, dt, U0, Ud0, Fex, varargin)
%HHTAMARCH conducts time series marching using HHTA
  
  opts = struct('alpha', 0, 'beta', 1/4, 'gamma', 1/2, ... % Newmark
      'reletol', 1e-6, 'etol', 1e-6, 'utol', 1e-6, 'rtol', 1e-6, ...
      'Display', 'progress', ...  % can be 'iter', 'progress', 'both', 'waitbar'
      'ITMAX', 10);
  if length(varargin)==1
      nflds = fieldnames(varargin{1});
      for i = 1:length(nflds)
        opts.(nflds{i}) = varargin{1}.(nflds{i});
      end
  end  
  
  a = opts.alpha;
  b = opts.beta;
  g = opts.gamma;
  
  
  U0  = reshape(U0, m.Ndofs, 1);
  Ud0 = reshape(Ud0, m.Ndofs, 1);
  
  Z1 = m.M + (1+a)*g*dt*m.C + (1+a)*b*dt^2*m.K;
  Z2 = m.M - (1+a)*(1-g)*dt*m.C - (1+a)*(0.5-b)*dt^2*m.K;
  Z3 = (1+a)*dt*m.K;
  
  T = T0:dt:T1;    Nt = length(T);
  U = zeros(m.Ndofs, Nt);  Ud = U;  Udd = U;
  U(:, 1) = U0;  Ud(:, 1) = Ud0;   
  
  % Initialize acceleration
  [Fnl, ~, ~, m] = m.NLFORCE(T0, U0, Ud0, T0-dt);
  Udd(:,1) = m.M\(Fex(T0)-m.C*Ud0-m.K*U0-Fnl);
  clear U0 Ud0
  
  if strcmp(opts.Display, 'waitbar')
      wb = waitbar(T(1)/T1, sprintf('Progress: %.e/%.e', T(1), T1), ...
          'createcancelbtn', "setappdata(gcbf, 'interrupt', true)");
  end
  for i=2:Nt
      %% Explicit Predictor
      Udd(:, i) = Udd(:, i-1);
      
      %% Corrector Iterations
      [FnlP, dFnldu, dFnldud, ~] = m.NLFORCE(T(i-1)+(1+a)*dt, ...
          U(:, i-1) + (1+a)*dt*Ud(:, i-1) + (1+a)*dt^2*((.5-b)*Udd(:, i-1)+b*Udd(:,i)), ...
          Ud(:, i-1) + (1+a)*dt^2*((1-g)*Udd(:, i-1)+g*Udd(:, i)), T(i-1));
      % Residual, Jacobian, and Update
      R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
          (FnlP-Fnl) - (Fex(T(i)+(1+a)*dt)-Fex(T(i-1)));  
      J = Z1 + (1+a)*(b*dt^2*dFnldu + g*dt*dFnldud);
      du = -J\R;
      % Error norms
      e = abs(R'*du);
      r = mean(R.^2);
      u = mean(du.^2);
      
      e0 = e;
      r0 = r;
      u0 = u;
      it = 0;
      
      flag = 8*(e/e0<opts.reletol) + 4*(e<opts.etol) + 2*(r<opts.rtol) + ...
          1*(u<opts.utol);
      if strcmp(opts.Display, 'iter') || strcmp(opts.Display, 'both')
          fprintf('ITN, E, E/E0, r, du\n%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, flag);
      end
      while (flag<=8) || (it==0)
          Udd(:, i) = Udd(:, i) + du;
          it = it+1;
          
          [FnlP, dFnldu, dFnldud, ~] = m.NLFORCE(T(i-1)+(1+a)*dt, ...
              U(:, i-1) + (1+a)*dt*Ud(:, i-1) + (1+a)*dt^2*((.5-b)*Udd(:, i-1)+b*Udd(:,i)), ...
              Ud(:, i-1) + (1+a)*dt^2*((1-g)*Udd(:, i-1)+g*Udd(:, i)), T(i-1));
          % Residual, Jacobian, and Updates
          R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
              (FnlP-Fnl) - (Fex(T(i)+(1+a)*dt)-Fex(T(i-1)));
          J = Z1 + (1+a)*(b*dt^2*dFnldu + g*dt*dFnldud);
          du = -J\R;
          % Error norms
          e = abs(R'*du);
          r = mean(R.^2);
          u = mean(du.^2);
          
          flag = 8*(e/e0<opts.reletol) + 4*(e<opts.etol) + 2*(r<opts.rtol) + ...
              1*(u<opts.utol);
          if strcmp(opts.Display, 'iter') || strcmp(opts.Display, 'both')
              fprintf('%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, flag);
          end
          
          if it>opts.ITMAX
              flag = 0;
              break;
          end
      end
      
      if flag == 0 || any(~isfinite(abs(U(:, i))))
          disp(sprintf('No Convergence/Non finite march at %f s : Returning', T(i)))
%          keyboard
          
          U = U(:, 1:i-1);
          Ud = Ud(:, 1:i-1);
          Udd = Udd(:, 1:i-1);
          T = T(:, 1:i-1);
          break;
      end
      
      %% Update States
      Ud(:, i) = Ud(:, i-1) + dt*((1-g)*Udd(:, i-1)+g*Udd(:, i));
      U(:, i) = U(:, i-1) + dt*Ud(:, i-1) + dt^2*((0.5-b)*Udd(:, i-1)+b*Udd(:, i));
      
      [Fnl, ~, ~, m] = m.NLFORCE(T(i), U(:, i), Ud(:, i), T(i-1));
      
      if strcmp(opts.Display, 'progress') || strcmp(opts.Display, 'both')
          fprintf('---------------------------------------------------\n');
          fprintf('%.4e/%.4e %.4e\n', T(i), T1, dt);
          fprintf('---------------------------------------------------\n');
      end
      
      %% Check for kill
      if strcmp(opts.Display, 'waitbar')
          waitbar(T(i)/T1, wb, sprintf('Progress: %e/%e', T(i), T1))
          
          if (~ishandle(wb))
            break;
          elseif getappdata(wb, 'interrupt')
            delete(wb);
       
            U   =   U(:, 1:i);
            Ud  =  Ud(:, 1:i);
            Udd = Udd(:, 1:i);
            T   =   T(:, 1:i);       
            return;
          end
      end
  end
  if strcmp(opts.Display, 'waitbar')  
    waitbar(1.0, wb, 'COMPLETED!');
    delete(wb);
  end
end
