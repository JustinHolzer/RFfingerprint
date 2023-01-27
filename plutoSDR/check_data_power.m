%% make size(X) = #samples_per_frame x 2 x 1 x #Frames
clear all
load('NNdata_2Plutos_10000_Mar2022.mat')
X = double(X);
[samps, b2, c1, numFrames] = size(X);

for i = 1: numFrames
    frame = complex(X(:,1,1,i),X(:,2,1,i));
    pwr(i) = frame'*frame;
end

plot(pwr)