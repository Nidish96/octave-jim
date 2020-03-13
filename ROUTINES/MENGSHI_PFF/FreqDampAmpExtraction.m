function [Freq,Damp,Amp,Time]=FreqDampAmpExtraction(data,t)
%%
% note: the calculated frequency and damping are only based on the upper amplitude

% input
% data: reverse filtered data
% t: time

% output
% Freq: raw instantaneous frequency (Hz)
% Damp: raw instantaneous damping
% Amp: raw instantaneous amplitude
% Time: time, of the same length as the above three vectors

%% Peak Finding and Fitting
[X_max,T_max,~,~,~,~]=PeakFinding(data,t);
%% Calculation of Frequency
[Freq_up]=FreqCal(T_max);
%% Calculation of Damping
[Damp,Freq,Amp,Time]=DampCal_c(Freq_up,X_max(1:end-1),T_max(1:end-1));

end

