% FFT and IFFT examples

% creating signal
clc;
clear; % clears all variables
close all % close all windows
F1=20;
F2=100;
Fs = 1000;   % Sampling frequency
Ts = 1/Fs;        % Sampling time
t=0:Ts:.3;
x=sin(2*pi*F1*t)+sin(2*pi*F2*t);
L = numel(x);
noisy_x = x+.5*randn(1,L);  % Signal Corrupted with Zero-Mean Random Noise

% FFT
NFFT = pow2(nextpow2(L));  % equivalent to 2^(nextpow2(L)) 
              % {nextpow2 is function that returns exp of next higher power
              % of 2}
% NFFT = 301;
fft_x = fft(x,NFFT);        % take FFT
fft_noisy_x = fft(noisy_x,NFFT);
p1 = fft_x.*conj(fft_x);    % mult by complex conjigut to get power (p)
p2 = fft_noisy_x.*conj(fft_noisy_x);
f = Fs/2*linspace(0,1,NFFT/2+1); % determine the vector of frequencies

% IFFT
ifft_x = ifft(fft_x,NFFT);
ifft_noisy_x = ifft(fft_noisy_x,NFFT);

% plotting
figure(1);
subplot(4,2,1);
plot(t,x,'.-');
axis([0 .3 -2.5 2.5])
title('Original signal: sin(2*pi*20*t)+sin(2*pi*100*t)')
xlabel('time(sec)');
ylabel('x(t)');
grid on

subplot(4,2,2);
plot(t,noisy_x,'.-');
axis([0 .3 -2.5 2.5])
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time(sec)');
ylabel('Noisy x(t)');
grid on

subplot(4,2,3);
plot(f,2*abs(fft_x(1:NFFT/2+1)),'r.-')
title('Single-Sided Amplitude Spectrum of signal')
xlabel('freq(Hz)');
grid on

subplot(4,2,4);
plot(f,2*abs(fft_noisy_x(1:NFFT/2+1)),'r.-')
title('Single-Sided Amplitude Spectrum of noisy signal')
xlabel('freq(Hz)');
grid on

subplot(4,2,5);
plot(f,4*p1(1:NFFT/2+1),'g.-')
title('Power spectral density of original signal')
xlabel('freq(Hz)');
grid on

subplot(4,2,6);
plot(f,4*p2(1:NFFT/2+1),'g.-')
title('Power spectral density of noisy signal')
xlabel('freq(Hz)');
grid on

subplot(4,2,7);
plot((0:1/(NFFT-1):1)*NFFT*Ts,real(ifft_x),'k.-')
axis([0 .512 -2.5 2.5])
title('IFFT')
xlabel('time(sec)');
grid on

subplot(4,2,8);
plot((0:1/(NFFT-1):1)*NFFT*Ts,real(ifft_noisy_x),'k.-')
axis([0 .512 -2.5 2.5])
title('IFFT')
xlabel('time(sec)');
grid on

