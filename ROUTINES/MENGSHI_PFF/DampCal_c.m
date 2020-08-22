function [zt,zt_fn,zt_A,zt_t]=DampCal_c(fn,A,t)
%% use central difference to calculate damping ratio
% input (same lengths)
% fn: frequency
% A: amplitude
% t: time

% output
% zt: damping ratio
% zt_fn, zt_A, zt_t: the corresponding frequency, amplitude, and time
%%
kc1=1;
for loop1=2:length(A)-1
    Ad(kc1)=(log(A(loop1+1))-log(A(loop1-1)))/...
        (t(loop1+1)-t(loop1-1));
    kc1=kc1+1;
end

kc2=2;
for loop2=1:length(Ad)
    zt(loop2,1)=-Ad(loop2)/(2*pi*fn(kc2));
    kc2=kc2+1;
end
zt_t=t(2:end-1);
zt_A=A(2:end-1);
zt_fn=fn(2:end-1);
end




