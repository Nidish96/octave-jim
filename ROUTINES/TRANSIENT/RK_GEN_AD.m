function [time, y, z] = RK_GEN_AD(func, trange, IC, ICz, pars)
%RK_GEN_AD integrates generic function with static evolving parameter
%of the form,
%   [y', z] = func(t, y, z)
% USAGE:
%   [time, y, z] = RK_GEN_AD(func, trange, IC, ICz, pars);
% INPUTS:
%   func    : function handle as above
%   trange  : 1x2 time range
%   IC      : Initial conditions in the dynamic states
%   ICs     : Initial conditions in the static states
%   pars    : structure with parameters
%       a   : Butcher table a
%       b   : Butcher table b
%       bs  : Butcher table bs
%       c   : Butcher table c
%    abstol : absolute tolerance between steps
%    pow    : power for step adaptation
%    maxstep: maximum time step size

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % RK 45 Butcher Tableau
% % pars.a = [0 0 0 0 0 0; 
% %           1/4 0 0 0 0 0; 
% %           3/32 9/32 0 0 0 0; 
% %           1932/2197 -7200/2197 7296/2197 0 0 0;
% %           439/216 -8 3680/513 -845/4104 0 0;
% %          -8/27 2 -3544/2565 1859/4104 -11/40 0];
% % pars.b = [16/135 0 6656/12825 28561/56430 -9/50 2/55];
% % pars.bs = [25/216 0 1408/2565 2197/4104 -1/5 0];
% % pars.c = [0 1/4 3/8 12/13 1 1/2];
% % % Step size controls
% % pars.abstol = 1e-6;
% % pars.pow = 1/4;
% % pars.maxstep = 1e-3;
% % % Display
% % pars.Display = 'on';
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


    a = pars.a;
    b = pars.b;
    bs = pars.bs;
    c = pars.c;
    abstol = pars.abstol;
    pow = pars.pow;
    maxstep = pars.maxstep;
    minstep = 1e-6;
    
    ord = length(c);
    dt = pars.maxstep;
    Np = fix((trange(end)-trange(1))/dt+1);
    if Np<101
        Np = 101;
        dt = (trange(end)-trange(1))/(Np-1);
    end
    time = linspace(trange(1), trange(end), Np)';
    y = zeros(Np, length(IC));
    y(1, :) = IC;
    z = zeros(Np, length(ICz));
    z(1, :) = ICz;
    
    k  = zeros(ord, length(IC));
    kz = zeros(ord, length(ICz));
    n  = 1;
    while time(n)<trange(end)
        k = zeros(size(k));
        for i=1:ord
            % Option 1: Do not update z between RK steps
            % [k(i,:), kz(i,:)] = func(time(n)+c(i)*dt, (y(n,:)+dt*a(i,:)*k)', z(n,:)');
            
            % Option 2: Update z consistently (Better)
            [k(i,:), kz(i,:)] = func(time(n)+c(i)*dt,...
                (y(n,:)+dt*a(i,:)*k)', (z(n,:)+a(i,:)*(kz-z(n,:)))');
        end
        en = dt*norm((b-bs)*k);
%         en = dt*norm((b-bs)*k(:, 1:size(k,2)/2));  % Only error in displacement update
        er = abstol/en;
        if en==0
            en = dt*eps;
            er = abstol/en;
            chi = 50;
        else
            chi = er^pow;
        end
        
        if chi>0.8
            y(n+1,:) = y(n,:) + dt*b*k;
            
            % z(n+1, :) = kz(end,:);  % Update z as per last step
            z(n+1,:) = z(n,:) + b*(kz-z(n,:));  % Update z consistently (Better)
            
            time(n+1) = time(n)+dt;
            n = n+1;
            
            if strcmp(pars.Display,'on')
                fprintf("t=%e dt=%e en=%e\n",time(n),dt,en);
            end
        end
        dt = min(maxstep,chi*dt);
    end
    time(n:end) = [];
    y(n:end,:) = [];
    z(n:end,:) = [];
end