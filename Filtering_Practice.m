% Filtering practice

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
