% Practice on FFT and IFFT

% creating signal
clc;
clear; % clears all variables
close all % close all windows

F1=input('Please input first frequency component(Hz):');
F2=input('Please input second frequency component(Hz):');
Fs = 1000;   % Sampling frequency
Ts = 1/Fs;        % Sampling time
L = input('Please input number of samples:');
t=(0:L-1)*Ts;
x=sin(2*pi*F1*t)+sin(2*pi*F2*t)+.5*randn(1,numel(t));  % Signal Corrupted with Zero-Mean Random Noise 

% FFT
NFFT = input('Please input number of fft points (NFFT):');
fft_x = fft(x,NFFT);
p1 = fft_x.*conj(fft_x);
f = Fs/2*linspace(0,1,NFFT/2+1);

% IFFT
NIFFT = input('Please input number of ifft points (NIFFT):');
ifft_x = ifft(fft_x,NIFFT);

% plotting
figure(1);
subplot(311);
plot(t,x,'.-');
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time(sec)');
ylabel('x(t)');
grid on

subplot(312);
plot(f,2*abs(fft_x(1:NFFT/2+1)),'r.-')
title('Single-Sided Amplitude Spectrum of signal')
xlabel('freq(Hz)');
ylabel('Amplitude');
grid on

subplot(313);
plot((0:L-1)*Ts,real(ifft_x(1:L)),'k.-')
title('IFFT')
xlabel('time(sec)');
ylabel('Amplitude');
grid on