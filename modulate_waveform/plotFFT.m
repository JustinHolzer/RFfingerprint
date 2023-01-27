%
%
function plotFFT(x, Fs, figNum,plotColor,flag_fullFFT,nFft)
%x - input signal
%Fs - sampling freq (Hz)
%figNum (optional) - figure number (allows you to call function multiple times and add
          %traces to given figNum
%plotColor (optional)- color of trace on plot (e.g. ='.-r';)
%flag_fullFFT (optional)- shift fft and plots full FFT spectrum plot(fftshift(x)) 
                %(if flag_fullFFT = 0; then only plots "positive freqs"


if(~exist('nFft'))
    nFft = 2^nextpow2(length(x));
end
xFft = fft(x ,nFft);
x_freq = linspace(0,Fs/2,nFft/2);

if(~exist('figNum'))
    figure
elseif(figNum==-1)
    figure
else
    figure(figNum);
end

if(~exist('plotColor'))
    plotColor = '.-';
end

if(~exist('flag_fullFFT'))
    flag_fullFFT = 1;
end

if(~flag_fullFFT)    
    %plot half (real part) of spectrum
    %plot(x_freq,abs((xFft(1:nFft/2))),plotColor); 
    semilogy(x_freq,abs((xFft(1:nFft/2))),plotColor); 
    xlabel('freq (Hz)');
else
    %plot full (real and imag) spectral content
    x2_freq = linspace(-Fs/2,Fs/2,nFft);
    %plot(x2_freq,abs(real(modSigFft)),plotColor);
    if(flag_fullFFT>1)
        x2_freq_noshift = linspace(0,Fs,nFft);
        %plot(x2_freq_noshift,abs(((xFft))),plotColor);
        semilogy(x2_freq_noshift,abs(((xFft))),plotColor);
    else
        %plot(x2_freq,abs((fftshift(xFft))),plotColor);
        semilogy(x2_freq,abs((fftshift(xFft))),plotColor);
    end
    xlabel('frequency (Hz)');
end

hold on;