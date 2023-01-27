%% loadData4MTM.m   also loads pluto data for NN
% repurposed from WISPR
% Generate input data
% Encode with Convolutional Encoder
% modulate with MSK (maybe GMSK) modulation
%
% references: https://www.mathworks.com/help/comm/ug/estimate-ber-for-hard-and-soft-decision-viterbi-decoding.html
%
clear all
numTrainingFramesPerTxr = 500;  %number of collects per file =  numCollects = 500;
frameLength_vec = [2500, 2500, 2778, 2500, 2500,2500,2500];
flag_data_vec = [0,1,2,3,4,5,0.5];
for i = 4
flag_data = flag_data_vec(i)
frameLength = frameLength_vec(i);
if(flag_data==0) %MSKdata
    rx_folder = 'NNdata_MSKdata/'
    rx_filename = [rx_folder, 'filenames_MSKdata26-Apr-2022.mat'];
elseif(flag_data==0.5) %MSKdata - but choose data randomly from frame
    rx_folder = 'NNdata_MSKdata/'
    rx_filename = [rx_folder, 'filenames_MSKdata26-Apr-2022.mat'];
    frameLenShort = 4000; %500
    numTrainingFramesPerTxr = 500;
    file_NNdata = [rx_folder, ['NN_MSK_8Txrs_4000x2x1x500_26Apr.mat']]
elseif(flag_data==1) %
    rx_folder = 'NNdata_MSKones/'
    rx_filename = [rx_folder, 'filenames_MSKones26-Apr-2022.mat'];
elseif(flag_data==2) %
    rx_folder = 'NNdata_5GUL/'
    rx_filename = [rx_folder, 'filenames_5GUL26-Apr-2022.mat'];
elseif(flag_data==3) %
    rx_folder = 'NNdata_5GDL/'
    rx_filename = [rx_folder, 'filenames_5GDL26-Apr-2022.mat'];
elseif(flag_data==4) %
    rx_folder = 'NNdata_5GDL_w/'
    rx_filename = [rx_folder, 'filenames_5GDL_w26-Apr-2022.mat'];
elseif(flag_data==5) %
    rx_folder = 'NNdata_5GDL_w/'
    rx_filename = [rx_folder, 'filenames_5GDL_w26-Apr-2022.mat'];
    frameLenShort = 20; %500
    numTrainingFramesPerTxr = 100;
end
load(rx_filename)
num_Txrs=length(filenames); % 8 Txrs  
X = []; y = [];
preamble = definePreamble();
numMissedFrames = 0;
for i = 1:length(filenames)
    file = filenames{i}; disp(file)
    load([rx_folder,file]);
    
    x_all = double(rxdata);

    plutoCollectSamps = 5556;  %see "frameLength2" in "rx_mod_GMSK_5G_pluto_revb.m"
    idx_start = 1;
    for k = 1:numTrainingFramesPerTxr
        x_this = x_all(idx_start:idx_start+plutoCollectSamps-1);
        % find start of preamble
        if(flag_data==0.5)
            x = chooseSignalRandom(x_this, frameLenShort);
            X = [X; x];
            y = [y; i]; %transmitter classification number
        else
            [x, num_lag] = findSignal(x_this, preamble);
            mtmOut = pmtm(x);
            %plotFFT(x,Fs); title('FFT of signal')
            if(length(mtmOut)==4096)
 X = [X; mtmOut];
 y = [y; i];
            else
                numMissedFrames = numMissedFrames + 1;
            end
%                 y = [y; i]; %transmitter classification number
%             if(length(x) ~= frameLength) %2639)% %missed Frame (didn't detect preambles properly). Each pluto collect saved this many samples
%                 disp(['wrong number of symbols found in data,  see line 54','length of x is: ',num2str(length(x))])
% 
%                 if(length(x) == 2639 & (flag_data==2))
%                     x = [x; zeros(139,1)];
%                     X = [X; x];
%                     y = [y; i];
%                 elseif(length(x) == 2639 & (flag_data==3))
%                     x = x(1:end-139);
%                     X = [X; x];
%                     y = [y; i];
%                 else
%                     numMissedFrames = numMissedFrames + 1;
%                 end
%             else%   did not miss a Frame
%                 if(exist('frameLenShort'))
%                     x = x(1:frameLenShort);
%                 end
%                 X = [X; x];
%                 y = [y; i]; %transmitter classification number
%             end
        end
    end
end

%% format data for Neural Network
X = X(:);
X = [real(X) imag(X)];
%% make size(X) = #samples_per_frame x 2 x 1 x #Frames
if(exist('frameLenShort'))
    frameLength = frameLenShort;
end
frameLengthMTM = 4096;
X = permute(reshape(X,[frameLengthMTM, numTrainingFramesPerTxr*num_Txrs-numMissedFrames, 2, 1]), [1 3 4 2]);


%% shuffle the data 
[~,~,~,numFrames] = size(X);
[~, idx] = sort(rand(numFrames,1));
Xshuffled = X(:,:,:,idx);
X = Xshuffled;
yShuffled = y(idx);
y = yShuffled;

totalFrames = num_Txrs*numTrainingFramesPerTxr-numMissedFrames;
numTrainingFrames = totalFrames*0.8;
numTrainingFramesPerTxr =  numTrainingFrames/num_Txrs;
numTestFrames = totalFrames*0.1;
numValFrames = totalFrames*0.1;
xTrainingFrames = X(:,:,:,1:numTrainingFrames);
xValFrames = X(:,:,:,(numTrainingFrames+1):(numTrainingFrames+numValFrames));
xTestFrames = X(:,:,:,(numTrainingFrames+1+numValFrames):(numTrainingFrames+numValFrames+numTestFrames));
yTrain = y(1:numTrainingFrames);
yVal = y((numTrainingFrames+1):(numTrainingFrames+numValFrames));
yTest = y((numTrainingFrames+1+numValFrames):(numTrainingFrames+numValFrames+numTestFrames));
if(~exist('file_NNdata'))
    file_NNdata = [rx_folder, ['MSK_rand_',num2str(num_Txrs),'Txrs_2500x2x1x500_26Apr.mat']]
end
frameLength = frameLengthMTM;
save(file_NNdata,'X','xTrainingFrames','xValFrames','xTestFrames','y','yTrain','yVal','yTest','frameLength','num_Txrs','numTrainingFramesPerTxr')
end
    %%%  TBD ....
    % the freq offset is just for my interest to understand what is
    % happening with SDR
%     % find coarse freq offset (could be worth knowing...)
%     [rxSig_coarse, cFreqOffset] = findCoarseFreqOffset(rxSig, Fs, preamble);
%     % estimate fine freq/phase offset.
%     [x_ff_corrected, freq_offset] = findFineFreqOffset(rxSig_coarse, length(preamble), Fs);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%      Functions       %%%%%%%%%%%%%%%%%%%%%%
function x = chooseSignalRandom(x_this, frameLen)
    maxStartIdx = length(x_this)-frameLen-1;
    startIdx = floor(rand(1)*maxStartIdx); %figure(34); plot(startIdx,1,'x'); hold on; 
    x = x_this(startIdx+1:startIdx+frameLen);
end
%%%%%%%%%%% function  %%%%%%%%%%%%%%%
function loadChirpToNNFrames(file_NNdata, filenames,numTrainingFramesPerTxr)
%% load data to X (from all files)
num_Txrs = length(filenames);
frameLength = 160;
training_frac = 0.8;
X = []; y = [];
for i = 1:num_Txrs
    load(filenames{i})
    %if    % check for enough data ...  
    %(numTrainingFramesPerTxr*frameLength)>length(rxdata);fprintf('error'):
    stpIdx = 0;
    for kk = 1 : numTrainingFramesPerTxr
        strtIdx = stpIdx + 1;
        stpIdx = strtIdx + frameLength - 1;
        dataFrame = rxdata(strtIdx:stpIdx);
        X = [X; dataFrame];
        y = [y; i];
        %fprintf('.');
    end
    fprintf(' completed i= %i',i);
end
end

function loadPlutoDataNN(file_NNdata, filenames,numTrainingFramesPerTxr)
%% load data to X (from all files)
num_Txrs = length(filenames);
frameLength = 160;
training_frac = 0.8;
X = []; y = [];
for i = 1:num_Txrs
    load(filenames{i})
    %if    % check for enough data ...  
    %(numTrainingFramesPerTxr*frameLength)>length(rxdata);fprintf('error'):
    stpIdx = 0;
    for kk = 1 : numTrainingFramesPerTxr
        strtIdx = stpIdx + 1;
        stpIdx = strtIdx + frameLength - 1;
        dataFrame = rxdata(strtIdx:stpIdx);
        X = [X; dataFrame];
        y = [y; i];
        %fprintf('.');
    end
    fprintf(' completed i= %i',i);
end
end
function [x_payload n_lag] = findSignal(x, preamble)
    % x - input signal.  findSignal assumes there are two instances of the
    % preamble.  findSignal returns the signal/payload data between the two
    % preambles
    [C lags] = xcorr(x,preamble);
    [Cmax, indx] = max(C);
    n_lag = lags(indx)+1;
    xNoLag = x(n_lag:end);
    
    %figure; plot(abs(C))
    %%  find other peak(s)...
    C2 = C;
    n = 139*3;
    C2(indx-n:indx+2*n) = zeros(3*n+1,1);
    [Cmax2, indx2] = max(C2);
    n_lag2 = lags(indx2)+1;
    
    if(n_lag2<n_lag)
        n_start = n_lag2+length(preamble);  %remove preamble
        n_stop = n_lag-1;
    elseif(n_lag2>n_lag)
        n_start = n_lag+length(preamble);
        n_stop = n_lag2-1;
    else
        error(['preambles not found - see function(findSignal)', num2str(dbstack)])
    end
    x_payload = x(n_start:n_stop);
    if(0)
        figure; %subplot(2,1,1)
        plot(real(preamble),'x-b'); hold on
        plot(real(x(n_lag:(n_lag+ length(preamble)-1))),'o-r')
    end
    if(0)
        figure; subplot(2,1,1)
        n=20;
        plot(real(preamble(1:n)),'x-b'); hold on
        subplot(2,1,2) %n=length(preamble);
        plot(real(x(n_lag:(n_lag+n-1 ))),'o-r')
    end
end
function preamble = definePreamble(numPreambleBits)
    %% define preamble
    if(~exist('numPreambleBits'))
        numPreambleBits = 139*2;
    end
    zcRoot = 25;
    preamble = zadoffChuSeq(zcRoot,numPreambleBits/2);
    preamble = [preamble; preamble];
end
function [xfreqFreqCorrected, coarseFreqOffset] = findCoarseFreqOffset(x, Fs, preamble)
    % assumes that a rough detection algorithm has been performed and Y
    % contains the reference waveform (preamble) preamble = refWaveform
    N = length(x);
    freq_search_range = 10e3;
    fshifts_hz = [-1*freq_search_range:2e3:freq_search_range]; %-50e3:5e3:50e3;  %% note: this should be wider than any potential oscillator mismatch
    F = zeros(N,length(fshifts_hz)); %help Matlab go faster
    % Pre-generate frequency corrections
    for k = 1:length(fshifts_hz)
        f = fshifts_hz(k); %[Hz] (potential) freq offset to be tested 
        t =(0:N-1)';
        F(:,k) = exp(1i*2*pi*f*t/Fs);
        Yf(:,k) = x.*F(:,k);
    
        [test,lag] = xcorr(Yf(:,k),preamble);
        tmp = test(find(lag==0):end);
        Yfcorr(:,k) = tmp;
        [maxVal(k), maxIndx(k)] = max(tmp);
        %figure(3); plot(abs(tmp),'x');
    end
    %%find max row from correlated Y with various freq offsets 
    %find which freq_shift gives highest correlation
    [tmp,ind2_max] = max(maxVal);
    coarseFreqOffset = fshifts_hz(ind2_max)
    xfreqFreqCorrected = Yf(:,ind2_max);
end

% fine frequency offset (reference: frequencyOffsetSSB_mmsesync function)
function [x_ff_corrected, freq_offset] = findFineFreqOffset(x, preamble_len, Fs)
    %inputs: x: preamble with repeated structure = [pream pream] where length(pream) = preamble_len/2
    %        Fs: sampling frequency
    if(0) %usage example
        f = 12345; %[Hz]
        x=preamble;
        t = (0:length(x)-1);
        x = x(:).';  % make sure x is row vector    
        x_freq_offset = x.*exp(1i*2*pi*f*t/Fs);
        freq_offset = findFineFreqOffset(x_freq_offset, length(preamble), Fs)
    end

    x = x(:).';  % make sure x is row vector    
    x1 = x(1:preamble_len/2); 
    x2 = x(preamble_len/2 + 1 : preamble_len); 
    x_dot = x1*x2';
    x_angle = angle(x_dot);
    freq_ratio = angle(x_dot)/(2*pi*length(x1));
    freq_offset = freq_ratio*Fs;  

    t =(0:length(x)-1);
    F_shift = exp(1i*2*pi*freq_offset*t/Fs);  %plotFFT(rxSig, Fs, 22)
    %rxSig_noShift = rxSig; %debug
    x_ff_corrected = x.*F_shift;                   %plotFFT(rxSig, Fs, 22)
end

