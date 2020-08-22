clc
clear all
addpath('../../ROUTINES/HARMONIC/')

rng(1)
fex = @(t) rand(size(t));

fsamp = 2^18;
dt = 1/fsamp;

t = 0:dt:1;
fscl = sqrt(length(t));
famp = 1;
ft = famp*fscl*fex(t);

% ft = wgn(1, length(t), 80);  % 40: 1

for i = 1:length(t)
    ft(i) = wgn(1, 1, 40);
end

[freqs, Ff] = FFTFUN(t', ft');

figure(1)
clf()
plot(t, ft)
xlabel('Time (s)')
ylabel('Force')

figure(2)
clf()
semilogy(freqs, abs(Ff), '.'); hold on
semilogy(freqs, famp*ones(size(freqs)), 'k--')

xlabel('Frequency (Hz)')
ylabel('Force')