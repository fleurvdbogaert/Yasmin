%% TM stage 2 - KHFAC Processing 
% Build Convolutional Neural Network 
% Y. (Yasmin) Ben Azouz 
% 4559843 

% Toolboxes required: 
% 1. Deep Learning Toolbox (dividerand function)
% 2. Fixed-point designer 
%% 1.   Load data 
load data2 
raw = data2(2,:) ; 
lab = data2(3,:) ; 

%% 2.   Normalize data using z-score
%data2(4,:) = cell(size(data2(1,:))) ; 
% for i = 1:length(raw)
%     raw{:,i} = normalize(raw{:,i},'zscore');
% end
%% 3.   Divide data into training, validation and test set (0.8 - 0.1 - 0.1)
[trainIdx,valIdx,testIdx] = dividerand(length(raw),0.8,0.1,0.1);

Xtrain = raw(trainIdx);
Xval = raw(valIdx);
Xtest = raw(testIdx);

Ytrain = lab(trainIdx);
Yval = lab(valIdx);
Ytest = lab(testIdx);
%% 4.   Divide each sequence into parts of equal length A, each fragment being 
%       shifted a preset number of samples from the previous fragment

% shift_train = 1;
% shift_test = 1;
% [Xequal_train, Yequal_train] = divide(Xtrain, Ytrain, A, shift_train);
% [Xequal_val, Yequal_val] = divide(Xval, Yval, A, shift_test);
% [Xequal_test, Yequal_test] = divide(Xtest, Ytest, A ,shift_test);

%% 5.  Concatenate data in validation and test set
% Xequal_val_together = [];
% for i=2:length(Xequal_val)
%   add2 = vertcat(Xequal_val{i});  
%   Xequal_val_together = [Xequal_val_together, add2];
% end
% 
% Yequal_val_together = [];
% for i=2:length(Yequal_val)
%   add5 = vertcat(Yequal_val{i});  
%   Yequal_val_together = [Yequal_val_together, add5];
% end
% 
% Xequal_test_together = [];
% for i=2:length(Xequal_test)
%   add3 = vertcat(Xequal_test{i});  
%   Xequal_test_together = [Xequal_test_together, add3];
% end
% 
% Yequal_test_together = [];
% for i=2:length(Yequal_test)
%   add6 = vertcat(Yequal_test{i});  
%   Yequal_test_together = [Yequal_test_together, add6];
% end

%%
numSequence = size(raw{1,1},1);
numHiddenUnits = 60;
InitialLearnRate = 0.0011113;
MaxEpochs = 2;
MiniBatchSize = 30;

layers = [...
    sequenceInputLayer(numSequence,'Normalization', 'zscore')
    bilstmLayer(numHiddenUnits,'OutputMode','sequence')
    fullyConnectedLayer(2)
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

