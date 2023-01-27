%%reference: https://www.mathworks.com/help/supportpkg/plutoradio/ug/repeated-waveform-transmitter.html
clear all;
close all;

%% Setup Pluto SDRs for single PC
TxNum_all =['00';'04'];
TxNum = TxNum_all(1,:);
RxNum = TxNum;
flag_plot = 1;
% RxDevice = sdrdev('Pluto');
% RxDevice.RadioID = 'usb:0';
% TxDevice = sdrdev('Pluto');
% TxDevice.RadioID = 'usb:1';

%rx = sdrrx('Pluto','RadioID','usb:0');
%tx = sdrtx('Pluto','RadioID','usb:1');

 numSamps = 10e3;
 Fs = 7680000;
 f = 1e6; % frequency for sine wave at baseband (note: f<Fs/2)
 t = [0:numSamps-1]/Fs;
 x = exp(1i*2*pi*f*t);
if(0)
    figure;  plot(real(x(1:1e2)),'x-')
    hold on; plot(imag(x(1:1e2)),'o-r')
end
    
    %% Transmitter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %tx_filename = '50_UL_FRC_1Frames_5MHzBW_16QAM.mat' % gives rx data with no drop-outs (drop-outs mean where abs(rxdata)~0)
%     %tx_filename = '4_DL_FRC_100Frames_5MHzBW_QPSK.mat'; % gives rx data with some drop-outs (where abs(rxdata)~0)
%     %% chirp files:
%     if(file == 1); tx_filename = 'chirp_25e4BW_10kSamps.mat';
%         BW = '250kHz'; samps = '10kSamps';
%     elseif(file == 2); tx_filename = 'chirp_25e5BW_10kSamps.mat';
%         BW = '2.5MHz'; samps = '10kSamps';
%     elseif(file == 3); tx_filename = 'chirp_25e4BW_20kSamps.mat';
%         BW = '250kHz'; samps = '20kSamps';
%     elseif(file == 4); tx_filename = 'chirp_25e5BW_20kSamps.mat';
%         BW = '2.5MHz'; samps = '20kSamps';
% %     elseif(file == 5); tx_filename = 'chirp_25e5BW_50kSamps.mat'; % 50k sample files appear to be longer than the capture period of the Pluto
% %         BW = '2.5MHz'; samps = '50kSamps'; numSamps = 50e3;
% %     elseif(file == 6); tx_filename = 'chirp_25e4BW_50kSamps.mat';
% %         BW = '250kHz'; samps = '50kSamps'; numSamps = 50e3;
%     elseif(file == 5); tx_filename = 'chirp_25e4BW_15kSamps.mat';
%         BW = '250e3'; samps = '15kSamps'; numSamps = 15e3;
%     elseif(file == 6); tx_filename = 'chirp_25e5BW_15kSamps.mat';
%         BW = '2.5e6'; samps = '15kSamps'; numSamps = 15e3;
%     
%     end
%     %tx_filename = [tx_filename(1:end-4),'2.mat'];
%     load(tx_filename);
%     
    waveform = x(:);
    cntrFreq = 3560e6; %tells what freq for Pluto to upconvert baseband signal to
    tx_gain = 0; %<=0
    % Create a transmitter System object for the AD936x-based radio hardware and set desired radio settings.
    tx = sdrtx('Pluto', 'CenterFrequency', cntrFreq, 'BasebandSampleRate', Fs, 'Gain',tx_gain);
    % Send waveform to the radio and repeatedly transmit it.
    txWaveform = int16(floor((2^15*waveform)+0.5));
    transmitRepeat(tx,txWaveform);

    %% Receiver  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(1)
        %filenames = {};
        numFiles = 4;
        for k = 1:numFiles
            %rx_filename = ['chirp_Tx',TxNum,'_Rx',RxNum,'_',BW,'_',samps,'_',num2str(k),'.mat'];
            rx_filename = ['sine_Tx',TxNum,'_Rx',RxNum,'_',num2str(k),'.mat'];
            %filenames = [filenames; rx_filename];
            rx_gain = 20; % -4 to 71 Radio receiver gain in dB. (reference: https://www.mathworks.com/help/supportpkg/plutoradio/ref/comm.sdrrxpluto-system-object.html)
            rx = sdrrx('Pluto');
            rx.BasebandSampleRate = Fs;
            rx.CenterFrequency = cntrFreq;
            %rx.Gain = rx_gain
            %rx.SamplesPerFrame = rx.BasebandSampleRate * 2;

            % transmitRepeat(tx, data_in)
            [rxdata,datavalid,overflow] = rx();
            numCollects = 100;
for ii = 1:numCollects
    % transmitRepeat(tx, data_in)
    [tmp_rxdata,datavalid,overflow] = rx();
    rxdata = [rxdata; tmp_rxdata];
    p_time = 0.015; pause(p_time);
    fprintf('.');
    %fprintf(['length(rxdata): ', num2str(length(rxdata)),'  pause(',num2str(p_time),')\n'])
end

            fprintf('Collecting samples...');
            p_time = 3; pause(p_time);
            fprintf('Done!\n');
            %fprintf(['length(rxdata): ', num2str(length(rxdata)),'  pause(',num2str(p_time),')\n'])

            release(rx);
            ts = 1/Fs*([1:length(rxdata)]-1);
            if(flag_plot)
%                 figure; %subplot(2,1,1)
%                 plot(ts,real(rxdata),'x-'); title(rx_filename,'Interpreter', 'none')
                %subplot(2,1,2); plot(ts,abs(double(rxdata)),'x-')
                plotFFT(rxdata,Fs); title(rx_filename,'Interpreter', 'none')
            end

            %save(rx_filename, "rxdata", "cntrFreq","Fs","datavalid","overflow","rx_gain","ts")
        end
    end

    release(tx)
if(0)
    filenames = {};
    TxNum_all =['00';'04'];
    numFiles = 4;
    for i = 1:length(TxNum_all)
        TxNum = TxNum_all(i,:);
        RxNum = TxNum;
        for k =1:numFiles
            rx_filename = ['sine_Tx',TxNum,'_Rx',RxNum,'_',num2str(k),'.mat'];
            filenames = [filenames; rx_filename];
        end
    end
    save(['filenames_sine',date,'.mat'],"filenames")
    %load(['filenames_sine',date,'.mat']);
end