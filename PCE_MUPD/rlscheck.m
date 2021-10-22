clc
clear all

%% 
emg = 0.2;
f = @(x) tanh((x-0.5)*10);
fn = @(x) tanh((x-0.5)*10)+(rand(size(x))-0.5)*emg;

Nx = 100;
x = linspace(0, 1, Nx+1)';  x(1) = [];
plot(x, f(x), '-', x, fn(x), '.')

%% Least-Squares Line
Phi = [x ones(Nx,1)];
Y = fn(x);

Theta_ls = (Phi'*Phi)\(Phi'*Y);

figure(1)
clf()
plot(x, Y, '.', x, Phi*Theta_ls, '-')

%% Recursive Least-Squares 
lamb = 0.85;
alph = 100;
theta_0 = Theta_ls;
% theta_0 = [0; -1];
P_0 = 1000*eye(length(theta_0));
Thetas = zeros(2, Nx+1);

% theta_0 = theta_0(1);
% P_0 = P_0(1,1);

figure(1)
clf()
for t=1:Nx    
%     P_t = P_0 - (P_0*Phi(t,:)'*Phi(t,:)*P_0)/(1+Phi(t,:)*P_0*Phi(t,:)');

	P_t = (P_0 - (P_0*Phi(t,:)'*Phi(t,:)*P_0)/(lamb+Phi(t,:)*P_0*Phi(t,:)'))/lamb;  % Exponential Forgetting (lambda)    
    theta_t = theta_0 + P_t*Phi(t,:)'*(Y(t)-Phi(t,:)*theta_0);
    
%     P_t = (P_0 - (P_0*x(t)*x(t)*P_0)/(lamb+x(t)*P_0*x(t)))/lamb;  % Exponential Forgetting (lambda)    
%     theta_t = theta_0 + P_t*x(t)*(Y(t)-x(t)*theta_0);
    
    figure(1)
    plot(x, Y, '.', x, Phi*Theta_ls, '-', ...
        x(t), Y(t), 'ko', x, Phi*theta_t, '-')
%     plot(x, Y, '.', x, Phi*Theta_ls, '-', ...
%         x(t), Y(t), 'ko', x, x*theta_t, '-')    
    ylim(1.5*[-1 1])
    pause(0.02)
    
    % Update
    theta_0 = theta_t;
    P_0 = P_t;
end