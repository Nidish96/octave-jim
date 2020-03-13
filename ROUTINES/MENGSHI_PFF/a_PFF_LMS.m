%% calculate freq, damp, and amp from LMS data
clear all;close all;clc;

%% Input (need to be changed properly)
str1='Run_1_14Mar_High';  % name of the data file
DOF=1; % the column to be analyzed 
fn=178.3; % central frequency in filtering (Hz)
n_bp=4; % order of the bandpass filter
bw_bp=60; % bandwidth of frequency in the bandpass filter
TimeDuration=1.5; % the time length to be analyzed 
% NOTE: TimeDuration depends on noise (not too long, can be changed according to the calculated results below)
%% Read Raw Data and Filter
load([str1])
[x,t_raw,dt,fs,f]=ReadLMS(Signal_0,Signal_1); % raw data of all DOFs
data_raw=x(:,DOF);
[xR,~,t] = RevForFilt(data_raw,t_raw,fn,n_bp,bw_bp);
[ X,X_abs,Fz1 ] = FFT_simp(data_raw,fs);
[ Y,Y_abs,Fz2 ] = FFT_simp(xR,fs);
%%
figure % FFT of raw data and filtered data
semilogy(Fz1,X_abs)
hold on;
semilogy(Fz2,Y_abs,'r')
xlim([0,3000])
legend('raw','filtered')
title('FFT of raw data and filtered data')
figure % time history of raw data and filtered data
subplot(2,1,1)
plot(t_raw,data_raw)
title('Raw data')
subplot(2,1,2)
plot(t,xR,'r')
title('Filtered data')
xlabel('time')
%% Cut Data Length
temp1=find(t<=TimeDuration);temp1=temp1(end);
xR=xR(1:temp1);t=t(1:temp1);

%% Calculation
[Freq,Damp,Amp,Time]=FreqDampAmpExtraction(xR,t);
%% Figure for Frequency (Backbone Curve)
figure
semilogx(Amp,Freq,'.');
xlabel('Amplitude');
ylabel('Frequency (Hz)');

%% Figure for Damping
figure
semilogx(Amp,Damp*100,'.');
xlabel('Amplitude');
ylabel('Damping Ratio (%)');

disp('Done!!!!!!!!!!!!!!!!!!!')



