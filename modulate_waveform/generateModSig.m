clear all; close all
flag_plot = 1;

Fs = 7680000;

for i = 1:6
    file=i;
    if(file == 1); tx_filename = 'chirp_25e4BW_10kSamps.mat';
        BW = 25e4; samps = '10kSamps'; numSamps = 10e3;
    elseif(file == 2); tx_filename = 'chirp_25e5BW_10kSamps.mat';
        BW = 2.5e6; samps = '10kSamps'; numSamps = 10e3;
    elseif(file == 3); tx_filename = 'chirp_25e4BW_20kSamps.mat';
        BW = 25e4; samps = '20kSamps'; numSamps = 20e3;
    elseif(file == 4); tx_filename = 'chirp_25e5BW_20kSamps.mat';
        BW = 2.5e6; samps = '20kSamps'; numSamps = 20e3;
    elseif(file == 5); tx_filename = 'chirp_25e4BW_15kSamps.mat';
        BW = 250e3; samps = '15kSamps'; numSamps = 15e3;
    elseif(file == 6); tx_filename = 'chirp_25e5BW_15kSamps.mat';
        BW = 2.5e6; samps = '15kSamps'; numSamps = 15e3;
    end
    t = [0:numSamps-1]/Fs;
    %f0 = 20e3;
    %x = chirp(t,f0,t(end),signalBW,'linear',0,'complex'); % same as doing hilbert transform on a real valued
    % signal, this complex-valued chirp shifts the FFT to only real-valued
    % frequencies.  alternatively  x = chirp(t,f0,t(end),signalBW); x=hilbert(x);
    x = chirp(t,-1*BW/2,t(end),BW/2,'linear',0,'complex');
    %txTime = numSamps/Fs;
    if(flag_plot)
        plotFFT(x,Fs); title('FFT of x')
        figure;
        pspectrum(x,Fs,'spectrogram')%,'TimeResolution',0.1, 'OverlapPercent',99,'Leakage',0.85)
    end

    %% define preamble
    numPreambleBits = 139*2;
    zcRoot = 25;
    preamble = zadoffChuSeq(zcRoot,numPreambleBits/2);
    preamble = [preamble; preamble];
    xpre = [preamble(:); x(:)]; %column vector with preamble and chirp
    if(flag_plot)
        plotFFT(xpre,Fs)
    end


    %% format data to waveStruct for Pluto SDR
    waveStruct.waveform = xpre;
    waveStruct.Fs = Fs;
    waveStruct.signalBW = BW;
    waveStruct.numSamps = numSamps;

    %tx_filename = [tx_filename(1:end-4),'2.mat'];
    save(tx_filename, 'waveStruct')

    tx_file2 = ['noPreamble_',tx_filename];
    waveStruct.waveform = x;
    waveStruct.preamble = preamble;
    save(tx_file2, 'waveStruct')

end

% if(0) %debug and experimentation...
%     f1 = 1;
%     f2 = signalBW;
%     x1 = exp(1i*pi*(f2-f1)*t.^2);
%     plotFFT(x1,Fs)
%
%     x3 = hilbert(x);
%
%     plotFFT(x,Fs); title('FFT of x')
%     plotFFT(x3,Fs); title('Hilbert Transform')
%
%     [a ind]=find(imag(x3)<-1.2);
%     figure
%     plot(t,imag(x3)); hold on
%     plot(t(ind),imag(x3(ind)),'x-r')
% end
