function [freq]=FreqCal(t)
%% use backward difference to calculate frequency
kc=1;
for loop1=2:length(t)
    freq(kc,1)=1/(t(loop1)-t(loop1-1));
    kc=kc+1;
end