%%reference: https://www.mathworks.com/help/supportpkg/plutoradio/ug/repeated-waveform-transmitter.html
clear all;
clc;
close all;


load('DL_FRC_100Frames_5MHzBW_QPSK.mat');
waveform = waveStruct.waveform;
sampleRate = waveStruct.Fs;
cntrFreq = 3560e6;
% Create a transmitter System object for the AD936x-based radio hardware and set desired radio settings.
tx = sdrtx('Pluto', 'CenterFrequency', cntrFreq, 'BasebandSampleRate', sampleRate);
% Send waveform to the radio and repeatedly transmit it.
txWaveform = int16(floor((2^15*waveform)+0.5));
transmitRepeat(tx,txWaveform);
%release(tx)
if(1)
    rx = sdrrx('Pluto');
    rx.BasebandSampleRate = sampleRate;%4e6;
    rx.CenterFrequency = cntrFreq;
    %rx.SamplesPerFrame = rx.BasebandSampleRate * 2;
    
   % transmitRepeat(tx, data_in)
    [rxdata,datavalid,overflow] = rx();
    
    fprintf('Collecting samples...');
    pause(5);
    fprintf('Done!\n');
    
    release(rx);
    plot(real(rxdata),'x-')
end

