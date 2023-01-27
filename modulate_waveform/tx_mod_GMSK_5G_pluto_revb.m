%% %% tx_mod_GMSK_5G_pluto_revb.m
%%reference: https://www.mathworks.com/help/supportpkg/plutoradio/ug/repeated-waveform-transmitter.html
clear all;
%clc;
close all;

flag_data = 0  %0,1,2,3;      0=MSK(rand()), 1=MSK(ones()) 2=5G Uplink data 3=5G DL data

if(flag_data == 0 | flag_data ==1)
    [waveform, Fs, numBits, nsamp] = makeMSKwaveform(flag_data);
    % debug waveform = [waveform; zeros(1500,1)]; % debug 
elseif(flag_data ==2)
    numSamps = 250*10;
    fn = '50_UL_FRC_1Frames_5MHzBW_16QAM.mat';
    [waveform, Fs] = load5Gwaveform(fn, numSamps);
elseif(flag_data ==3)
    %fn = '1_DL_FRC_1Frames_5MHzBW_64QAM.mat';
    fn = '2_DL_FRC_10Frames_5MHzBW_QPSK.mat';
    numSamps = 250*10; %payload data size (not including preamble) +139*2;
    [waveform, Fs] = load5Gwaveform(fn, numSamps);
end

%% Transmitter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cntrFreq = 2415e6; %3560e6;
tx_gain = 0; %<=0
% Create a transmitter System object for the AD936x-based radio hardware and set desired radio settings.
tx = sdrtx('Pluto', 'CenterFrequency', cntrFreq, 'BasebandSampleRate', Fs, 'Gain',tx_gain,'RadioID','usb:0');
connectedRadios = findPlutoRadio;
tmp = connectedRadios.SerialNum %pulls out first serial number (which is used by this script for the transmitter
if(tmp == '104473222a87000702002c001f4dd520b7') %serial number of Rxr (which should not be selected)
    release(tx)
    usb_num = 'usb:1';
    disp('release and reconnect, was connected tx to usb:0 serial ...d520b7') %, tmp, 'which is usb:1 ...  reconnecting to usb:1'])
    tx = sdrtx('Pluto', 'CenterFrequency', cntrFreq, 'BasebandSampleRate', Fs, 'Gain',tx_gain,'RadioID',usb_num);
end
% Send waveform to the radio and repeatedly transmit it.
txWaveform = int16(floor((2^15*waveform)+0.5));
transmitRepeat(tx,txWaveform);
connectedRadios = findPlutoRadio
connectedRadios.RadioID
connectedRadios.SerialNum

if(0)
    release(tx)
end


%% %%  Functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [waveform, Fs, numBits, nsamp] = makeMSKwaveform(flag_data)
    %BW = 500e6; %[Hz]
    %symRate = BW;
    numBits = 250;
    nsamp = 10; %MSK samples per symbol
    Fs = 7680000;
    
    rng(123); %set random number generator seed
    
    %% generate input data (payload)
    if(flag_data==0)
        data = randi([0 1],numBits,1);
    elseif(flag_data ==1)
        data = ones(numBits,1);
    end
    
    %% define preamble
    numPreambleBits = 139*2;
    zcRoot = 25;
    preamble = zadoffChuSeq(zcRoot,numPreambleBits/2);
    preamble = [preamble; preamble];
    
    %% MSK Modulate
    dataMSK = mskmod(data,nsamp);  %note: length(dataMSK) = 8*length(data_before_MSKModulation)
    % add preamble after modulation
    waveform = [preamble; dataMSK];
    
    %% 5G Modulation
    if(0)
        plotFFT(waveform,Fs); title('FFT of MSK data with preamble')
        plotFFT(dataMSK,Fs); title('FFT of MSK-modulated data')
        if(1); t = [0:numBits*nsamp-1]/Fs; figure; plot(t,real(dataMSK),'x-'); end
    end  
end
 z mnm                                          
function [waveform, Fs] = load5Gwaveform(fn, numSamps)
    load(fn);
    waveform5G = waveStruct.waveform;
    %figure; n=40; plot(phase(waveform5G(1:n)),'x-'); hold on; plot(phase(waveform5Gnorm(1:n)),'o-b')
    if(exist('numSamps'))
        n=10000;
        waveform5G = waveform5G(n+1:n+numSamps);
        %waveform5G = waveform5G(numSamps:2*numSamps-1);
    end
    waveform5G = 1.3*waveform5G/max(abs(waveform5G)); %normalize waveform
    %% define preamble
    numPreambleBits = 139*2;
    zcRoot = 25;
    preamble = zadoffChuSeq(zcRoot,numPreambleBits/2);
    preamble = [preamble; preamble];
    % add preamble after 5G signal
    waveform = [preamble; waveform5G];
    Fs = waveStruct.Fs;
    if(1)  %plots - debug and sanity checks
        plotFFT(waveform,Fs); title('FFT of 5G data with preamble')
        plotFFT(waveform5G,Fs); title('FFT of 5G data')
        if(1); t = [0:length(waveform5G)-1]; figure; plot(t,real(waveform5G),'x-'); title('real(waveform5G)'); end
        if(1); t = [0:length(waveform)-1]; figure; subplot(2,1,1); plot(t,real(waveform),'x-');title('real(waveform)'); subplot(2,1,2); plot(t,imag(waveform));title('imag(waveform)');end
        if(0); t = [0:length(preamble)-1]/Fs; figure; plot(t,real(preamble),'x-'); title('real(preamble)'); end
    end
end

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
    plot(x_freq,abs((xFft(1:nFft/2))),plotColor);
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
end
