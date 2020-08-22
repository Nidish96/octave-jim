%% forward filtered data and reversed filtered data
function [xR,xF,t,temp] = RevForFilt(xt,t,fn,n_bp,bw_bp)
%%
% input:
% xt: orignal signal
% t: time
% fn: natural frequency of the signal (Hz)
% n_bp: order of the bandpass
% bw_bp: bandwidth of the bandpass

% output:
% xR: reverse filtered signal
% xF: forward filtered signal
% t: return the time t>=0 (since experimental data may include data before t=0)
%%
fs = 1/(t(2)-t(1));
xt_rev = xt(end:-1:1);
fc_bp = [fn-bw_bp/2,fn+bw_bp/2 ]/fs*2;
ftype = 'bandpass';
[b,a] = butter(n_bp,fc_bp,ftype);
xF = filter(b,a,xt);
xR = filter(b,a,xt_rev);
xR = xR(end:-1:1);
temp = find(t>=0,1);
xF = xF(temp:end);
xR = xR(temp:end);
t = t(temp:end);
end