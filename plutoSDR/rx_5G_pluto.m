%%reference: https://www.mathworks.com/help/supportpkg/plutoradio/ug/repeated-waveform-transmitter.html
clear all;
clc;
%close all;

tx_filename = '50_UL_FRC_1Frames_5MHzBW_16QAM.mat'
load(tx_filename); 
sampleRate = waveStruct.Fs;
cntrFreq = 3560e6;

%% Receiver  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rx_filename = 'txPluto1_rxPluto0_txfile50.mat'
rx_gain = 10; % -4 to 71 Radio receiver gain in dB. (reference: https://www.mathworks.com/help/supportpkg/plutoradio/ref/comm.sdrrxpluto-system-object.html)
rx = sdrrx('Pluto');
rx.BasebandSampleRate = sampleRate;
rx.CenterFrequency = cntrFreq
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

save(rx_filename, "rxdata", "cntrFreq","sampleRate","datavalid","overflow") %,"rx_gain")


