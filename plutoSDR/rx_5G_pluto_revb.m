%% rx_5G_pluto_revb.m
%%reference: https://www.mathworks.com/help/supportpkg/plutoradio/ug/repeated-waveform-transmitter.html
clear all;
%clc;
%close all;
flag_tx = 0;
flag_rx = 1;

%% Setup 2 Pluto SDRs for single PC

% RxDevice = sdrdev('Pluto');
% RxDevice.RadioID = 'usb:0';
% TxDevice = sdrdev('Pluto');
% TxDevice.RadioID = 'usb:1';

%rx = sdrrx('Pluto','RadioID','usb:0'); 
%tx = sdrtx('Pluto','RadioID','usb:1');


%% Transmitter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tx_filename = '50_UL_FRC_1Frames_5MHzBW_16QAM.mat' % gives rx data with no drop-outs (drop-outs mean where abs(rxdata)~0)
%tx_filename = '4_DL_FRC_100Frames_5MHzBW_QPSK.mat'; % gives rx data with some drop-outs (where abs(rxdata)~0)
load(tx_filename);

waveform = waveStruct.waveform;
sampleRate = waveStruct.Fs;
cntrFreq = 3560e6;
 
if(flag_tx)
    tx_gain = 0; %<=0
    % Create a transmitter System object for the AD936x-based radio hardware and set desired radio settings.
    tx = sdrtx('Pluto', 'CenterFrequency', cntrFreq, 'BasebandSampleRate', sampleRate, 'Gain',tx_gain);
    % Send waveform to the radio and repeatedly transmit it.
    txWaveform = int16(floor((2^15*waveform)+0.5));
    transmitRepeat(tx,txWaveform);
    %release(tx)
end

%% Receiver  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flag_rx)
    rx_filename = 't4r0_txf50_10M_b.mat'  %'t2r2_txf50_10M.mat' 'test.mat' %
    rx_folder = 'pluto_data_1M/';
    rx_folder = 'C:\Users\JU25855\OneDrive - MIT Lincoln Laboratory\1.WPI\Fingerprinting_Research\5G_Toolbox_Waveforms\data_29Jan2022'; %C:\Users\JU25855\OneDrive - MIT Lincoln Laboratory\CS534AI\project\5G_Toolbox_Waveforms\data_29Jan2022\'
    rx_folderFile = [rx_folder, rx_filename];
    rx_gain = 0; % -4 to 71 Radio receiver gain in dB. (reference: https://www.mathworks.com/help/supportpkg/plutoradio/ref/comm.sdrrxpluto-system-object.html)
    rx = sdrrx('Pluto');
    rx.BasebandSampleRate = sampleRate;
    rx.CenterFrequency = cntrFreq;
    %rx.Gain = rx_gain 
    %rx.SamplesPerFrame = rx.BasebandSampleRate * 2;
    
    fprintf('Getting set to collect samples...\n');
    numCollects = 500;  %was 500
    rxdata = [];
    for ii = 1:numCollects
        % transmitRepeat(tx, data_in)
        [tmp_rxdata,datavalid,overflow] = rx();
        rxdata = [rxdata; tmp_rxdata];
        p_time = 0.015; pause(p_time);
        fprintf('.');
        %fprintf(['length(rxdata): ', num2str(length(rxdata)),'  pause(',num2str(p_time),')\n'])
    end
    fprintf('Done!\n');
    release(rx);
    ts = 1/sampleRate*([1:length(rxdata)]-1);
    figure; 
    plot(ts,real(rxdata),'x-')

    connectedRadios = findPlutoRadio
    connectedRadios.RadioID
    connectedRadios.SerialNum
    disp('runs on 29 Jan 2022 used Pluto0 as Receiver: Serial #: 104473222a87000702002c001f4dd520b7')
    save(rx_folderFile, "rxdata", "cntrFreq","sampleRate","datavalid","overflow","rx_gain","ts")
end

