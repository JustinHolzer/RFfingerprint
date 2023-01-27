%% rx_mod_GMSK_5G_pluto_revb.m
%%reference: https://www.mathworks.com/help/supportpkg/plutoradio/ug/repeated-waveform-transmitter.html
clear all;
%clc;
%close all;

flag_data = 4;  %0,1,2,3,4;      0=MSK(rand()), 1=MSK(ones()) 2=5G Uplink data, 3=5G DL data, 4= 5G DL over wireless channel
chooseTx = 1;; %1 through 8

flag_firstTime_createRxFilenames = 1; %=1 first time through; =0 every other time

TxNum_all =['14';'16';'15';'01';'02';'03';'10';'04'];  %4/21/22
TxNum = TxNum_all(chooseTx,:)
RxNum = '00';
Fs = 7680000;
cntrFreq = 2415e6; %3560e6;
frameLength2 = 2778*2;
%% Setup 2 Pluto SDRs for single PC
%% Transmitter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Receiver  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %rx_filename = 't4r0_txf50_10M_b.mat'  %'t2r2_txf50_10M.mat' 'test.mat' %
    %rx_folder = 'pluto_data_1M/';
    if(flag_data==0) %MSKdata
        rx_folder = 'NNdata_MSKdata/'        
        rx_filename = ['MSKdata_Tx',TxNum,'_Rx',RxNum,'_',date,'.mat'];
        %frameLength2 = 2778*2;
    elseif(flag_data==1) %
        rx_folder = 'NNdata_MSKones/'
        rx_filename = ['MSKones_Tx',TxNum,'_Rx',RxNum,'_',date,'.mat'];
        %frameLength2 = 2778*2;
    elseif(flag_data==2) %
        rx_folder = 'NNdata_5GUL/'
        rx_filename = ['5GUL_Tx',TxNum,'_Rx',RxNum,'_',date,'.mat'];
    elseif(flag_data==3) %
        rx_folder = 'NNdata_5GDL/'
        rx_filename = ['5GDL_Tx',TxNum,'_Rx',RxNum,'_',date,'.mat'];
    elseif(flag_data==4) %
        rx_folder = 'NNdata_5GDL_w2/'
        rx_filename = ['5GDL_Tx_w',TxNum,'_Rx',RxNum,'_',date,'.mat'];
    end
    if(flag_firstTime_createRxFilenames); mkdir(rx_folder); end
    rx_folderFile = [rx_folder, rx_filename];
    rx_gain = 0; % -4 to 71 Radio receiver gain in dB. (reference: https://www.mathworks.com/help/supportpkg/plutoradio/ref/comm.sdrrxpluto-system-object.html)
    rx = sdrrx('Pluto');
    rx.BasebandSampleRate = Fs; %sample Rate
    rx.CenterFrequency = cntrFreq;
    %rx.Gain = rx_gain 
    %rx.SamplesPerFrame = rx.BasebandSampleRate * 2;
    
    fprintf('Getting set to collect samples...\n');
    numCollects = 50;  %was 150, 500
    rxdata = [];
    for ii = 1:numCollects
        % transmitRepeat(tx, data_in)
        [tmp_rxdata,datavalid,overflow] = rx(); 

        rxdata = [rxdata; tmp_rxdata(1:frameLength2)];
        %p_time = 0.015; pause(p_time);
        fprintf('.');
        %fprintf(['length(rxdata): ', num2str(length(rxdata)),'  pause(',num2str(p_time),')\n'])
    end
    %disp(['number of data samples collected in one rx() = ',num2str(length(tmp_rxdata))])
    fprintf('Done!\n');
    release(rx);
    ts = 1/Fs*([1:length(rxdata)]-1);
            %%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOTS
            figure;   plot(ts*Fs,real(rxdata),'x-')
            plotFFT(rxdata(1:floor(length(rxdata)/10)),Fs)

    connectedRadios = findPlutoRadio
    connectedRadios.RadioID
    connectedRadios.SerialNum
    disp('runs on 21 Apr 2022 used Pluto0 as Receiver: Serial #: 104473222a87000702002c001f4dd520b7')
    save(rx_folderFile, "rxdata", "cntrFreq","Fs","datavalid","overflow","rx_gain","ts")  %% save data file 

    %% file with list of data file names
    if(flag_data==0) %MSKdata
        %rx_folder = 'NNdata_MSKdata/'        
        rxFilenames = ['filenames_MSKdata',date,'.mat'];
    elseif(flag_data==1) %
        %rx_folder = 'NNdata_MSKones/'
        rxFilenames = ['filenames_MSKones',date,'.mat'];
    elseif(flag_data==2) %
        %rx_folder = 'NNdata_5GUL/'
        rxFilenames = ['filenames_5GUL',date,'.mat'];
    elseif(flag_data==3) %
        %rx_folder = 'NNdata_5GDL/'
        rxFilenames = ['filenames_5GDL',date,'.mat'];
    elseif(flag_data==4) %
        %rx_folder = 'NNdata_5GDL_w/'
        rxFilenames = ['filenames_5GDL_w',date,'.mat'];
    end
    rxFilenames = [rx_folder,rxFilenames];
    if(flag_firstTime_createRxFilenames)  %initialize rxFile to store names of files
        filenames = {};
        save(rxFilenames,"filenames")
    end
    load(rxFilenames)
    filenames = [filenames; rx_filename]
    save(rxFilenames,"filenames")  
    %filenames = filenames(1:end-1) %% to erase last filename


