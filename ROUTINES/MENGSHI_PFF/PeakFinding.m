function [X_max,T_max,X_min,T_min,ind_max,ind_min]=PeakFinding(data,t)
%% Peak Finding and Fitting
% local peaks:  x_max;  local valleys:  x_min;  and their corresponding time :  t_max and t_min
[x_max,t_max,ind_max]=LocalPeak(data,t);
[x_min,t_min,ind_min]=LocalValley(data,t);
%% fitted peaks:  X_max;  corresponding time :  T_max
if data(1)==x_max(1)
    t_max_1=t_max(2:end);
    [X_max,T_max]=FitPV(data(2:end),t(2:end),t_max_1);
    T_max=[t_max(1);T_max];
    X_max=[x_max(1);X_max];
else if data(1)~=x_max(1)
        [X_max,T_max]=FitPV(data,t,t_max);
    end
end
%% fitted valleys:  X_min;  corresponding time :  T_min
if data(1)==x_min(1)
    t_min_1=t_min(2:end);
    [X_min,T_min]=FitPV(data(2:end),t(2:end),t_min_1);
    T_min=[t_min(1);T_min];
    X_min=[x_min(1);X_min];
else if data(1)~=x_min(1)
        [X_min,T_min]=FitPV(data,t,t_min);
    end
end

end

