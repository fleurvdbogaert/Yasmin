%% General Information
% Y. (Yasmin) Ben Azouz 
% 4559843
% July 2022 

%% Variables
DIR_IMPORT = [];        % when desired, set a default directory to select files from, otherwise use = [];
DIR_EXPORT = [];        % when desired, set a default directory the OUTPUT directory user interface will open up, otherwise use = [];
mkdir 'Figures'
mkdir 'Calculations'
FS_POT = 60000;        % [Hz] sampling frequency
FS_PRESS = 60000;

%% Input values
% C = 5111;            % Calibration factor
% C = 5420;            % calibration factor 150420
% C = 5859.4;          % Calibration factor 040320
C_press1 = 0;          % Calibration factor 16022022 = 28725;
C_press2 = 0;

TX  = 4;                % Amount of seconds of which the declining pressure slope upon inhibition will be calculated (seconds after start inhibition)
DS_FACT = 2;            % Downsampling factor
F_FILT = 18;            % [Hz] low-pass filter frequency for pressure 

% Set the number of stimulation blocks for low-frequency and high-frequency separately 
HF_BLOCKS = 1;
LF_BLOCKS = 0;

% Determine which pressure channels you want to evaluate
PRESS1 = 0;
PRESS2 = 1;
             
%% Execute
[SData, DIR_EXPORT] = load_filesV3(DIR_IMPORT, DIR_EXPORT);

fs_new = FS_POT /  DS_FACT; 

potential1 = SData.potential1;
pressure2 = SData.pressure2; 
%%

ROOT_DIR = uigetdir ; 
d = dir(fullfile(ROOT_DIR)) ; % and return all those matching files
d(1:3) = [] ;                 % remove '...', '...' and '.DS.files'
data = cell(5,numel(d)) ;     % create 3x35 cell array, top names bottom data  
for i=1:numel(d)
    data(1,i) = {d(i).name} ; 
    data(2,i) = {importdata(fullfile(ROOT_DIR,d(i).name)).stimulation(:,1)} ; %take the second (anal) pressure 

    absol = abs(data{2,i}) ; 

    m = movstd(absol,10000) ;
    data(3,i) = {m} ; 
    data(4,i) = {max(m)-min(m)} ; 
end 
%% plot 
k = 14 ; 

figure
subplot(2,1,1) ; plot(x,data{2,k}) ; title('raw signal')
subplot(2,1,2) ; plot(x,data{3,k}) ; title ('movstd')
%% 
[~,locs] = findpeaks(data{3,k},'MinPeakProminence',30); % 30 is goede threshold voor NoStim
%%
% make sure to base of of ONE simulation 
% removing section where no measurement was done
if isempty(locs) == 0 
    % start of measurement not correct 
    [zks,~] = find(~data{3,k}(1:locs(1))) ;
        if issorted(zks) == 1 
            data(5,k) = {data{3,k}(locs(1):end)} ;
        end
    % end of measurement not correct 
    [fks,~] = find(~data{3,k}(locs(end):end)) ;
        if issorted(fks,'monotonic') == 1 
            data(5,k) = {data{3,k}(1:locs(end))} ;
        end 
    % re-calculate peaks 
    [~,locs] = findpeaks(data{5,k},'MinPeakProminence',30);
    % re-calculate range 
    data(4,k) = {max(data{5,k})-min(data{5,k})} ;
elseif isempty(locs) == 1 
    % differentiate between stim and nostim through range 
    if data{4,k} < 10 % boundary for nostim range 
        % no stimulation given determined
        % What does this mean for the rest of the calculations? 
        % Do we leave the stim start and stop empty? 
    elseif data{4,k} > 10 % boundary for full stim? 

    end 
end 
%%
xz = linspace(1,length(data{5,k}),length(data{5,k})) ; 
figure 
plot(xz,data{5,k},...
    xz(locs),data{5,k}(locs),'rx') ; 
%%
figure(1); 
xa = linspace(1, length(a), length(a)) ; 
x = linspace(1, length(potential1), length(potential1)) ;
subplot(3,1,3) ; plot(xa,a) ; title('differential smooth')
subplot(3,1,1) ; plot(x,potential1) ; title('unfiltered')
subplot(3,1,2) ; plot(x,m) ; title ('integral')

%% plot findchangepts
TF = findchangepts(T,'Statistic','mean','MinThreshold',1000000) ;
 
%%
figure 
plot(x,T,...
    x(TF),T(TF),'rx') ; 

    % net_filter = designfilt('bandstopiir','FilterOrder',2, ...                       % create notch filter for 50 Hz (48-52 Hz)
%                            'HalfPowerFrequency1',49.9,'HalfPowerFrequency2',50.1, ...
%                            'DesignMethod','butter','SampleRate',fs_new);
% %fvtool(net_filter, 'FS', fs_new)                                       % Visualize the filter       
% potential1_filtered = filtfilt(net_filter, potential1);
%%
xb = linspace(1, length(b), length(b)) ; 
plot(xb,b)
