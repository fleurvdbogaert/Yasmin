%% TM stage 2 - KHFAC Processing 
% Build Convolutional Neural Network 
% Y. (Yasmin) Ben Azouz 
% 4559843 

% Toolboxes required: 
% 1. Deep Learning Toolbox (dividerand function)
% 2. Fixed-point designer 
%% 1.   Load data 
load data2 
data = data2 ; 
raw = data(2,:) ; 
lab = data(3,:) ; 

%% 2.   Normalize data using z-score
data(4,:) = cell(size(data2(1,:))) ; 
nor = raw ; 

A = 250;
for i = 1:length(raw)
    nor{:,i} = normalize(raw{:,i},'zscore');

    L = length(nor{:,i})/A;
    Ldes = A*ceil(L);
        for ii = length(nor{:,i}):Ldes
            nor{:,i}(ii) = 0;
        end
end
%% downsample 
raw2 = raw ; 
lab2 = lab ; 
for i = 1:size(raw,2) 
    raw2(:,i) = {downsample(raw{:,i},1000)} ; 
    lab2(:,i) = {downsample(lab{:,i},1000)} ; 
end 

%% Categorical names 
lab3 = lab2 ; 
valueset = [0 1] ; 
catnames = {'n/a' 'contraction'} ; 

for i = 1:size(lab2,2)
     B = categorical(lab2{:,i},valueset,catnames) ; 
     lab3(:,i) = {B} ; 
end 
%% 3.   Divide data into training, validation and test set (0.8 - 0.1 - 0.1)
[trainIdx,valIdx,testIdx] = dividerand(length(nor),0.8,0.1,0.1);

Xtrain = raw2(trainIdx)';
Xval = raw2(valIdx)';
Xtest = raw2(testIdx)';

Ytrain = lab3(trainIdx)';
Yval = lab3(valIdx)';
Ytest = lab3(testIdx)';
%%
Ytrain = cellfun(@(C) reshape(C, 1, []), Ytrain, 'UniformOutput',false);
%% 4.   Divide each sequence into parts of equal length A, each fragment being 
%       shifted a preset number of samples from the previous fragment

clear Xequal_train; clear Yequal_train; clear Xequal_test; clear Yequal_test;
shift_train = 1;
shift_test = 1;
[Xequal_train, Yequal_train] = divide(Xtrain, Ytrain, A, shift_train);
[Xequal_val, Yequal_val] = divide(Xval, Yval, A, shift_test);
[Xequal_test, Yequal_test] = divide(Xtest, Ytest, A ,shift_test);
%% 5. train
numSequence = size(raw2{1,1},1);
numHiddenUnits = 60;
InitialLearnRate = 0.0011113;
MaxEpochs = 2;
MiniBatchSize = 30;

layers = [...
    sequenceInputLayer(numSequence,'Normalization', 'zscore')
    bilstmLayer(numHiddenUnits,'OutputMode','sequence')
    fullyConnectedLayer(3)
    softmaxLayer
    classificationLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',MaxEpochs, ...
    'MiniBatchSize',MiniBatchSize, ...
    'InitialLearnRate',InitialLearnRate, ...
    'LearnRateDropPeriod',3, ...
    'LearnRateSchedule','piecewise', ...
    'GradientThreshold',1, ...
    'Plots','training-progress',...
    'shuffle','every-epoch',...
    'SequenceLength',500,...
    'Verbose',0,...
    'ValidationData',{Xval,Yval},...
    'ValidationFrequency',50,...
    'ValidationPatience',inf,...
    'DispatchInBackground',true );

net = trainNetwork(Xtrain, Ytrain, layers, options);
%% Onthouden 

layer = convolution2dLayer(filterSize,numFilters,Name,Value) ; % padding mogelijkheid 

