function [x_min,t_min,ind_min]=LocalValley(data,t)
%% find the local minimum of the data, and the corresponding time
dt=t(2)-t(1);
if data(1)<data(2) & data(1)<0
    kc=1;
    x_min(kc,1)=data(1);
    t_min(kc,1)=t(1);
    ind_min(kc,1)=1;
    kc=kc+1;
else
    kc=1;
end
for loop1=2:length(data)-1
    if data(loop1)<0 & data(loop1)<data(loop1+1) & data(loop1)<data(loop1-1)
        x_min(kc,1)=data(loop1);
        t_min(kc,1)=t(loop1);
        ind_min(kc,1)=loop1;
        kc=kc+1;
    end
end
if x_min(1)>x_min(2)
    x_min=x_min(2:end);
    t_min=t_min(2:end);
    ind_min=ind_min(2:end);
end

