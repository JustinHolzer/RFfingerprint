%%reference: https://www.mathworks.com/help/supportpkg/plutoradio/ug/repeated-waveform-transmitter.html
clear all;
clc;
%close all;
flag_2Plutos = 0;
if(flag_2Plutos)
    %% Setup Pluto SDRs for single PC
    RxDevice = sdrdev('Pluto');
    RxDevice.RadioID = 'usb:0';
    TxDevice = sdrdev('Pluto');
    TxDevice.RadioID = 'usb:1';
elseif(flag_2Plutos == 0)
    RxDevice = sdrdev('Pluto');
    RxDevice.RadioID = 'usb:0';
    RxDevice = sdrdev('Pluto');
    RxDevice.RadioID = 'usb:0';
end
%rx = sdrrx('Pluto','RadioID','usb:0'); 
%tx = sdrtx('Pluto','RadioID','usb:1');


%% Transmitter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tx_filename = '50_UL_FRC_1Frames_5MHzBW_16QAM.mat' % gives rx data with no drop-outs (drop-outs mean where abs(rxdata)~0)
%tx_filename = '4_DL_FRC_100Frames_5MHzBW_QPSK.mat'; % gives rx data with some drop-outs (where abs(rxdata)~0)
load(tx_filename); 

waveform = waveStruct.waveform;
sampleRate = waveStruct.Fs;
cntrFreq = 3560e6;
tx_gain = 0; %<=0
% Create a transmitter System object for the AD936x-based radio hardware and set desired radio settings.
tx = sdrtx('Pluto','RadioID','usb:1', 'CenterFrequency', cntrFreq, 'BasebandSampleRate', sampleRate, 'Gain',tx_gain);
% Send waveform to the radio and repeatedly transmit it.
txWaveform = int16(floor((2^15*waveform)+0.5));
transmitRepeat(tx,txWaveform);
%release(tx)

%% Receiver  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
    rx_filename = 't11_r10_txf50_2.mat'
    rx_gain = 10; % -4 to 71 Radio receiver gain in dB. (reference: https://www.mathworks.com/help/supportpkg/plutoradio/ref/comm.sdrrxpluto-system-object.html)
    rx = sdrrx('Pluto','RadioID','usb:0');
    rx.BasebandSampleRate = sampleRate;
    rx.CenterFrequency = cntrFreq;
    %rx.Gain = rx_gain 
    %rx.SamplesPerFrame = rx.BasebandSampleRate * 2;
    
   % transmitRepeat(tx, data_in)
    [rxdata,datavalid,overflow] = rx();
    
    fprintf('Collecting samples...');
    p_time = 4; pause(p_time);
    fprintf('Done!\n');
    %fprintf(['length(rxdata): ', num2str(length(rxdata)),'  pause(',num2str(p_time),')\n'])
    
    release(rx);
    ts = 1/sampleRate*([1:length(rxdata)]-1);
    figure; 
    plot(ts,real(rxdata),'x-')

    save(rx_filename, "rxdata", "cntrFreq","sampleRate","datavalid","overflow","rx_gain","ts")
end

