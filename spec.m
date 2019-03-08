function spec(x, Fs, nFFT)
% plots the spectra of input signal x
    fft_x = fft(x,nFFT);
    fft_x = fft_x(1:nFFT/2+1);
    
    f = Fs/2*linspace(0,1,nFFT/2+1);

    figure
    plot(f,fft_x.*conj(fft_x))
    xlabel('frequency'); ylabel('power')
end

