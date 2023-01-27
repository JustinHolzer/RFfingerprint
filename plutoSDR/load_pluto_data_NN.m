 %script: load_pluto_data_NN.m
% Purpose: this script is used to read and format the data collected from Pluto SDRs 
% so it can be processed by the Matlab Neural Network (see NeuralNetworkFingerprint.m)
clear all
numTrainingFramesPerTxr = 10000; %200; %62500; %125; %5000;
% filenames = {'txPluto1_rxPluto0_txfile50_1.mat';'txPluto2_rxPluto0_txfile50_1.mat';
%         'txPluto3_rxPluto0_txfile50_1.mat';'txPluto0_rxPluto0_txfile50.mat';
%         't10_r10_txf50_0.mat';'t14_r14_txf50_0.mat';'t15_r15_txf50_0.mat';
%         't16_r16_txf50_0.mat';'t17_r17_txf50_0.mat';'t18_r18_txf50_0.mat'};
addpath('pluto_data_1M')
%addpath 'C:\Users\JU25855\OneDrive - MIT Lincoln Laboratory\CS534AI\project\5G_Toolbox_Waveforms\pluto_data_1M'
addpath('C:\Users\JU25855\OneDrive - MIT Lincoln Laboratory\1.WPI\Fingerprinting_Research\5G_Toolbox_Waveforms\data_29Jan2022\data_mat_files')
%filenames = {'t0r0_txf50_10M.mat','t1r1_txf50_10M.mat','t2r2_txf50_10M.mat','t3r3_txf50_10M.mat'};
%filenames = {'t0r0_txf50_10M.mat','t1r1_txf50_10M.mat','t2r2_txf50_10M.mat','t3r3_txf50_10M.mat'};
filenames = findFilenames();
%num_Txrs = length(filenames);
dir_name = 'apr8_pluto_data/';
num_Txrs = 2;
mkdir(dir_name);
for i = 1:length(filenames)
    my_files{1} = filenames{i};
    for k = i+1:length(filenames)
        file_NNdata = [dir_name,'NNdata_',num2str(num_Txrs),'Plutos_',num2str(numTrainingFramesPerTxr),'_Apr2022_',num2str(i),num2str(k),'.mat']
        my_files{2} = filenames{k};
        loadPlutoDataNN(file_NNdata, my_files,numTrainingFramesPerTxr)
    end
end
%%%%%%%%%%% function  %%%%%%%%%%%%%%%
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


%% format data for Neural Network
X = X(:);
X = [real(X) imag(X)];
%% make size(X) = #samples_per_frame x 2 x 1 x #Frames
X = permute(reshape(X,[frameLength, numTrainingFramesPerTxr*num_Txrs, 2, 1]), [1 3 4 2]);


%% shuffle the data 
[~,~,~,numFrames] = size(X);
[~, idx] = sort(rand(numFrames,1));
Xshuffled = X(:,:,:,idx);
X = Xshuffled;
yShuffled = y(idx);
y = yShuffled;

totalFrames = num_Txrs*numTrainingFramesPerTxr;
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
save(file_NNdata,'X','xTrainingFrames','xValFrames','xTestFrames','y','yTrain','yVal','yTest','frameLength','num_Txrs','numTrainingFramesPerTxr')
end
%%%%%%%%%%% function  %%%%%%%%%%%%%%%
function files = findFilenames(folder, file_indx)
    if(~exist('folder'))
        folder = 'C:\Users\JU25855\OneDrive - MIT Lincoln Laboratory\1.WPI\Fingerprinting_Research\5G_Toolbox_Waveforms\data_29Jan2022\data_mat_files';
    end
    listing = dir(folder);
    if(~exist('file_nums'))
        file_indx = 1:length(listing);
    end
    %t0r0_txf50_10M.mat
    %maxNumFiles = 2; %length(listing);
    numFiles = 0;
    for i = 1:length(listing)
        %if((numFiles>=maxNumFiles) | (strcmp(listing(i).name,'.') | strcmp(listing(i).name,'..')))
        if( (strcmp(listing(i).name,'.') | strcmp(listing(i).name,'..')))
        else
            numFiles = numFiles + 1;
            files{numFiles} = listing(i).name;
        end
    end
end
