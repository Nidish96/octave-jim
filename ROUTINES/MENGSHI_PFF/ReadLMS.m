function [x,t,dt,fs,f]=ReadLMS(Signal_0,Signal_1)
% read time history response measured using LMS
% input:
% Signal_0, Signal_1: acceleration response of the measurement point
% and force response of the excitation point
% (the sequence does not matter)

% output:
% x: acceleration (time history response)
% t: time
% f: force (time history response)
    if Signal_1.y_values.quantity.label~='g'
        temp=Signal_1;
        Signal_1=Signal_0;
        Signal_0=temp;
    end
    x=Signal_1.y_values.values*9.8;
    dt=Signal_1.x_values.increment;
    fs=1/dt;
    t0=Signal_1.x_values.start_value;
    N=Signal_1.x_values.number_of_values;
    t=t0:dt:((N-1)*dt+t0);
    t=t';
    f=Signal_0.y_values.values;
end