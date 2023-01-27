function NeuralNetworkFingerprint(filename)
%%NeuralNetworkFingerprint.m
% This script reads in data that has been previously synthesized in Matlab 
% or collected from Software Defined Radios.  It is assumed the data is 
% formatted to fit in the Deep Neural Network (DNN) defined here.
% The data is then processed through an 18-layer DNN, including training,
% validating and testing the data. 


%%load x
if(~exist('filename'))
    %filename = ['X_dataMatrix_c','.mat']
    %filename = 'NNdata_10Plutos_125.mat';  %not enough training data
    %filename = 'NNdata_4Plutos_6250.mat'
    %filename = 'NNdata_4Plutos_62500.mat'
    %filename = 'NNdata_4Plutos_62500.mat'
    %filename = 'NNdata_9Plutos_200_Mar2022.mat'
    %filename = 'NNdata_9Plutos_5000_Mar2022.mat'
    filename = 'NNdata_2Plutos_10000_Mar2022.mat'
    addpath('data_simulation\')
    s = pwd;
    addpath([s,'\5G_Toolbox_Waveforms\'])
    addpath([s,'C:\Users\JU25855\OneDrive - MIT Lincoln Laboratory\1.WPI\Fingerprinting_Research\5G_Toolbox_Waveforms\'])
    %addpath('5G_Toolbox_Waveforms/pluto_data_1M')
end
load(filename);

%% Train the Neural Network
% This example uses a neural network (NN) architecture that consists of two convolutional and three fully connected layers. The intuition behind this design is that the first layer will learn features independently in I and Q. Note that the filter sizes are 1x7. Then, the next layer will use a filter size of 2x7 that will extract features combining I and Q together. Finally, the last three fully connected layers will behave as a classifier using the extracted features in the previous layers [1].
%reference: DesignRFFingerprintingDeepLearningNetworkExample.mlx


poolSize = [2 1];
strideSize = [2 1];
layers = [
  imageInputLayer([frameLength 2 1], 'Normalization', 'none', 'Name', 'Input Layer')
  
  convolution2dLayer([7 1], 50, 'Padding', [1 0], 'Name', 'CNN1')
  batchNormalizationLayer('Name', 'BN1')
  leakyReluLayer('Name', 'LeakyReLu1')
  maxPooling2dLayer(poolSize, 'Stride', strideSize, 'Name', 'MaxPool1')
  
  convolution2dLayer([7 2], 50, 'Padding', [1 0], 'Name', 'CNN2')
  batchNormalizationLayer('Name', 'BN2')
  leakyReluLayer('Name', 'LeakyReLu2')
  maxPooling2dLayer(poolSize, 'Stride', strideSize, 'Name', 'MaxPool2')
  
  fullyConnectedLayer(256, 'Name', 'FC1')
  leakyReluLayer('Name', 'LeakyReLu3')
  dropoutLayer(0.5, 'Name', 'DropOut1')
  
  fullyConnectedLayer(80, 'Name', 'FC2')
  leakyReluLayer('Name', 'LeakyReLu4')
  dropoutLayer(0.5, 'Name', 'DropOut2')
  
  fullyConnectedLayer(num_Txrs, 'Name', 'FC3')
  softmaxLayer('Name', 'SoftMax')
  classificationLayer('Name', 'Output')
  ]

%%  Configure the training options to use the ADAM optimizer with a mini-batch size of 256. 
%By default, 'ExecutionEnvironment' is set to 'auto', which uses a GPU for training if one is available. 
%Otherwise, trainNetwork uses a CPU for training. To explicitly set the execution environment, 
%set 'ExecutionEnvironment' to one of 'cpu', 'gpu', 'multi-gpu', or 'parallel'. 
%Choosing 'cpu' may result in a very long training duration.
trainNow = 1
if trainNow
  miniBatchSize = 256; %#ok<UNRCH>
  % Training options
  yVal = categorical(yVal);
  options = trainingOptions('adam', ...
    'MaxEpochs',100, ...
    'ValidationData',{xValFrames, yVal}, ...
    'ValidationFrequency',floor(numTrainingFramesPerTxr*num_Txrs/miniBatchSize/3), ...
    'Verbose',true, ...
    'L2Regularization', 0.0001, ...
    'InitialLearnRate', 0.0001, ...
    'MiniBatchSize', miniBatchSize, ...
    'ValidationPatience', 3, ...
    'Plots','training-progress', ...
    'Shuffle','every-epoch');
  
  % Train the network
  yTrain = categorical(yTrain);
  simNet = trainNetwork(xTrainingFrames, yTrain, layers, options);
  save('trainedNN.mat');  %save SimNet for future reference
else
  % Load trained network (simNet), testing dataset (xTestFrames and
  % yTest) and the used MACAddresses (generatedMACAddresses)
  myDir = pwd;
  addpath(myDir)
  load('rfFingerprintingSimulatedDataTrainedNN.mat',...
    'generatedMACAddresses',...
    'simNet',...
    'xTestFrames',...
    'yTest')
end

% Classify test frames and calculate the final accuracy off the neural network. 
% Obtain predicted classes for xTestFrames
yTestPred = classify(simNet,xTestFrames);

% Calculate test accuracy
yTest = categorical(yTest);
testAccuracy = mean(yTest == yTestPred);
disp("Test accuracy: " + testAccuracy*100 + "%")

% find confusion matrix for the test frames. As mentioned before, perfect classification accuracy is achieved with the synthetic dataset.
figure
cm = confusionchart(yTest, yTestPred);
cm.Title = 'Confusion Matrix for Test Data';
cm.RowSummary = 'row-normalized';

end
