

%% Sampling theorem:


    % pick a sampling frequency
    Fs = 120;

    % pick times of sampling
    t=0:1/Fs:1;
    w=12; % pick a signal frequency
    y = sin(w*t*2*pi);
    
    figure; plot(t, y, 'LineWidth',2); hold on
    
    % downsample the signal to 60 Hz
    t_60Hz = t(1:2:end);
    y_60Hz = y(1:2:end);
    
    plot(t_60Hz, y_60Hz, 'ro');
    
    plot(t_60Hz, y_60Hz, ':r', 'LineWidth',2);

    % downsample the signal to sub Nyquist
    t_10Hz = t(1:12:end);
    y_10Hz = y(1:12:end);
    
    plot(t_10Hz, y_10Hz, 'gd');
    
    plot(t_10Hz, y_10Hz, '-g', 'LineWidth',2);
    
    % they run ex1
    % they run ex1_practice
    
% =========================================================================    
%% +++++++++++++++++++++++++++++++++ FFT ++++++++++++++++++++++++++++++++++ 
% =========================================================================   

% creating signal
clc;
clear; % clears all variables
close all % close all windows

F1=20;
F2=100;
Fs = 1000;   % Sampling frequency
Ts = 1/Fs;        % Sampling time
t=0:Ts:.3;  % time vector

% Create signal
x=sin(2*pi*F1*t)+sin(2*pi*F2*t); % 2*pi to turn into radians

L = numel(x); % Length (number of elements)

noisy_x = x + 0.5*randn(1,L);  % Signal Corrupted with Zero-Mean Random Noise

% -------------------------------------------------------------------------
% FFT

NFFT = pow2(nextpow2(L));  % equivalent to 2^(nextpow2(L))
% NFFT = 301;

fft_x = fft(x,NFFT);        % take FFT
fft_noisy_x = fft(noisy_x,NFFT);
% FFT gives complex number

% let's look at the FFT:
fft_x(2)        % let's look at the 2nd element: 
                % it is a complex number - indicated by 'i'

plot(abs(fft_x)) % we can plot just the magnitude/amplitude with abs
% in plot, discard 2nd half cuz mirror image of the second half
hold on
plot(abs(fft_x).^2) % we can plot the magnitude squarred - power
% dot in the .^2 is to square each element
plot(fft_x.*conj(fft_x),'g+') % can plot the product with the complex conjigut
                              % are the same
hold off

p1 = fft_x.*conj(fft_x);    % mult by complex conjigut to get power (p)
p2 = fft_noisy_x.*conj(fft_noisy_x);

% f = Fs/2*linspace(0,1,NFFT/2+1); % determine the vector of frequencies
f = Fs/2*linspace(0,1,NFFT/2+1); %sampling rate

figure; plot(f,abs(fft_x(1:NFFT/2+1)))

% -------------------------------------------------------------------------

% IFFT (inverse FFT)
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

% Trailing 0s is due to number of points/components (want more components
% than is actually in your signal)
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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% play with ex2_FFT_IFFTPractice
%
% this will generate 1 second of signal
% 10, 50 hz, and enough input samples in 1 sec to satisy Nyquist (100) 
%  - pick number of FFT to be larger than num points (512)
%  - pick number of iFFT points to be the same (512)
% -> perfect reconstruction in the inverse FFT
%
% What if number of input samples is less than 100?
% What if number of FFT is smaller than number of samples?
% What if number of iFFT is smaller than FFT?
   
% MAKE SPEC function for plotting spectra

% lets make a signal
   Fs = 1000;
   t=0:1/Fs:1;
   x = sin(20*2*pi*t); plot(x)
   spec(x,Fs,512)
   
% what about a small square blip - say touched the apparatus and gave a
% static discharge.
    xs = zeros(1,numel(t));
    xs(round(numel(t)/2)-3:round(numel(t)/2)+3) = 40*ones(1,6+1);
    x = xs+ sin(10*2*pi*t); 
    plot(x)
    spec(x,Fs,512)  
    % why is this?
    % that does this mean for analysis of your data signals (have to get
    % rid of these)
    %   - how can you get rid of them: filtering, exclude traces
   
 % lets add some noise  
   x = sin(10*2*pi*t)+rand(1,numel(t)); plot(x)
   spec(x,Fs,512) % where is the low frequency part comming from? Rand has mean of 0.5
   
% what about a linear drift?
    x = sin(10*2*pi*t) + t - 0.5; plot(x)
    spec(x,Fs,512)
    % this has edge effects
    
% how can we get rid of these types of edge artifacts?
%  subtract mean and use a window that tapers
    x = sin(10*2*pi*t) + t - 0.5; plot(x)
   win = hanning(numel(x))';
    plot(x); hold on
    plot(win, 'r')
    plot(win.*x,'g')     % note that it still has some drift
    spec(win.*x,Fs,512)

    specCell({x,win.*x},Fs,512) % can see that the power of main component 
                        %is squashed, but edge effects are even more so.
 
 
% =========================================================================    
%% +++++++++++++++++++++++++++ Filtering ++++++++++++++++++++++++++++++++++ 
% =========================================================================  

Fs = 1000;
t =0:1/Fs:2;

% slow linear drift?
x = sin(10*2*pi*t) + t - 0.5; plot(x)
spec(x,Fs,512)
% this has edge effects (= first and last signal is not at center)

% how can we get rid of these types of edge artifacts?
%  subtract mean and use a window that tapers
win = hanning(numel(x))';
plot(x); 
hold on
plot(win, 'r')
plot(win.*x,'g')     % note that it still has some drift

spec(win.*x,Fs,512)

specCell({x,win.*x},Fs,512) % can see that the power of main component 
                    %is squashed, but edge effects are even more so.

% -------------------------------------------------------------------------

% low pass filter
%   what do we use for a low pass filter?
%   box-car or another kernel with smooth edges

Fs = 100; %sampling rate
t =0:1/Fs:2;
s1 = rand(1,numel(t))+.3*sin(2*pi*3*t); % make a signal w/randomness
plot(t,s1)

n = 5;
kern = ones(1,n)./n;        % box carr window (filter kernel)
                              % (want kernel to sum to 1)
s2 = conv(s1,kern,'same');  % convolve (multiply then sum); 
                              % same = same number of data points out as put in
figure; plot(t,s1,'b')
hold on
plot(t,s2,'r')

figure
subplot(2,1,1)
plot(t,s1); hold on
plot(t,s2,'r')
specCell({s1-mean(s1),s2-mean(s2)},Fs,64)

close 
clear

% -------------------------------------------------------------------------

% high-pass filter
%   what does a high-pass filter look like?

x = 1:200;
s1 = rand(1,200)-0.5;     % make a 0-mean random signal
s1_slow = s1 + sin(x/25); % add a slow component

n = 19;                  % kernal size
kern = -ones(1,n)./n;    % make a negative box-car
kern(10) = 1;            % set the middle element
s2 = conv(s1_slow,kern,'same'); % convole it

plot(s1_slow, 'b')
hold on
plot(s2,'r')
hold off


figure
subplot(2,1,1)
plot(s1); hold on
plot(s2,'r')
specCell({s1-mean(s1),s2-mean(s2)},Fs,64)

figure
plot(s1_slow); hold on   % plot the random signal with slow component
plot(s2,'r')             % plot the filtered trace - what happened?
plot(s1,'g')             % plot original trace
% is the filtered trace exactly the same as the original s1? (zoom in)

close 
clear

% -------------------------------------------------------------------------

% Mexican Hat
%  filters can do all sorts of things
%  Lets see a sum of 2 gaussians.

s1 = zeros(1,300);          % make a staircase signal
s1(100:200) = ones(1,101);
s1(201:300) = ones(1,100).*2;

figure
subplot(2,1,1)
plot(s1); hold on
hat = gausswin(100,10).*2 - gausswin(100,5).*0.8; % sum of 2 gaussians
plot(hat,'r') % this is our kernel
set(gca,'YLim',[-0.5,2.5])

subplot(2,1,2)
s2 = conv(s1, hat, 'same'); % convolution again
plot(s2)
% edge detector - contrast enhancement in vision



% phase lag
%
% all the kernels we have used so far have been symmetric
%  - what does that require? Hint, think about causal/non-causal
%  - knowledge in the future and the past
%  - you can do this after you have the signal
%  - * but what about the digitizer on your recording rig?
%    - that is analogue, so no memory

close 
clear

% -------------------------------------------------------------------------

% let's do the low-pass filter again
s1 = rand(1,200)+sin([1:200]/10);  % make a signal
plot(s1)

n = 10;
kern = ones(1,n)./n;        % box carr window
kern_NC = [0,0,kern,0,0];   % put 0 padding in symmetrically (not centered)
kern_C = [kern,0,0,0,0];    % put in 0 padding for future values (centered)
s2_NC = conv(s1,kern_NC,'same');  % convolve
s2_C = conv(s1,kern_C,'same');  % convolve

figure
% subplot(2,1,1)
plot(s1); hold on
plot(s2_NC,'r')
plot(s2_C,'g')
% you can see that there is a phase shift

% you have to be very careful about introducing phase shifts when
% filtering your data !!!!!




% how you build a filter in matlab.


 
 
function spec(x, Fs, nFFT)
% plots the spectra of input signal x
    fft_x = fft(x,nFFT);
    fft_x = fft_x(1:nFFT/2+1);
    
    f = Fs/2*linspace(0,1,nFFT/2+1); % just a funciton of #:inc:#

    figure
    plot(f,fft_x.*conj(fft_x))
    xlabel('frequency'); ylabel('power')
end


function specCell(xc, Fs, nFFT)
% plots the spectra of input cell array xc containing signals
figure
col_vec = {'b','g','r','m'};
    for i=1:numel(xc)
        x = xc{i};
        fft_x = fft(x,nFFT);
        fft_x = fft_x(1:nFFT/2+1);

        f = Fs/2*linspace(0,1,nFFT/2+1);

        plot(f,fft_x.*conj(fft_x), col_vec{mod(i-1,4)+1}); hold on
        xlabel('frequency'); ylabel('power')
    end
end

