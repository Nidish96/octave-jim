function [X_m,T_m]=FitPV(data,t,t_max_1)
%% fit the local extreme points based on the sampled ones
for loop1=1:length(t_max_1)
    temp(loop1)=find(t_max_1(loop1)==t); % index of local maxima
    x_temp(:,loop1)=data(temp(loop1)-1:temp(loop1)+1);
    t_temp(:,loop1)=t(temp(loop1)-1:temp(loop1)+1);
end
%%
for loop3=1:size(t_temp,2)
    A=[t_temp(1,loop3)^2 t_temp(1,loop3) 1;...
        t_temp(2,loop3)^2 t_temp(2,loop3) 1;...
        t_temp(3,loop3)^2 t_temp(3,loop3) 1;];
    y=[x_temp(1,loop3); x_temp(2,loop3); x_temp(3,loop3)];
    temp=A\y;
    a=temp(1);
    b=temp(2);
    c=temp(3);
    T_m(loop3,1)=-b/(2*a);
    X_m(loop3,1)=(4*a*c-b^2)/(4*a);    
end
