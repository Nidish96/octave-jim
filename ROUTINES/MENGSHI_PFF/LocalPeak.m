function [x_max,t_max,ind_max]=LocalPeak(data,t)
%% find the local maximum of the data, and the corresponding time
dt=t(2)-t(1);
if data(1)>data(2) & data(1)>0
    kc=1;
    x_max(kc,1)=data(1);
    t_max(kc,1)=t(1);
    ind_max(kc,1)=1;
    kc=kc+1;
else
    kc=1;
end
for loop1=2:length(data)-1
    if data(loop1)>0 & data(loop1)>data(loop1+1) & data(loop1)>data(loop1-1)
        x_max(kc,1)=data(loop1);
        t_max(kc,1)=t(loop1);
        ind_max(kc,1)=loop1;
        kc=kc+1;
    end
end
if x_max(1)<x_max(2)
    x_max=x_max(2:end);
    t_max=t_max(2:end);
    ind_max=ind_max(2:end);
end
