function [stim,pres,dir_outp] = loadModify(f_filt,fs_p,fs_s)
% Data is loaded into cells and modified. 
% %   Y. (Yasmin) Ben Azouz
%   NetID: 4559843
%   July 2022 
%   Only one file can be selected at a time 

%% Variables
% fs_new = fs_p / ds_fact;
au = 10 ;   % Conversion factor AU to cmH2O


%% Select input file and load data 
[name,path] = uigetfile('*.*', ...
    'Select ONE FILE to put in') ;
name = {name} ; 

% Display selected file in command window
if isequal(name{1,1},0)
   error('File selection canceled') ;
else 
   disp('File selected: '), disp(name{1,1}) ; 
   data = importdata(fullfile(path,name{1,1})) ;
end
%% Select output directory 
dir_outp = uigetdir('', 'Select OUTPUT directory');

%% Modify data 
% Stimulation 
num_stim = size(data.stimulation,2) ; % amount of stimulation measurements
stim = cell(2,num_stim) ;
for n = 1:num_stim
    stim(1,n) = {data.stimulation(:,n)} ; 
    stim(2,n) = stim(1,n) ; %no adjustments to data 
    stim(3,n) = {fs_p} ; 

    % Ask what to display
    answer = questdlg(strcat('Do you want to display stimulation channel ' ...
        , num2str(n),' ?'), ...
        'Select stimulation channels', ...
        'Yes','No','No');
    switch answer
        case 'Yes'
            disp(strcat('stimulation channel ', num2str(n),' will be displayed.'))
        case 'No'
            stim(:,n) = [] ; 
            disp(strcat('stimulation channel ', num2str(n),' will NOT be displayed.'))
    end   
end 
       
% pressure 
num_pres = size(data.pressure,2) ; % amount of pressure measurements 
pres = cell(2,num_pres) ; % first raw, second modified 
for nn = 1:num_pres
    pres(1,nn) = {data.pressure(:,nn)} ; 
    pres(3,nn) = {fs_s} ; 

    % Notch filter 50Hz
    net_filter = designfilt('bandstopiir','FilterOrder',2, ...                       % create notch filter for 50 Hz (48-52 Hz)
                       'HalfPowerFrequency1',49.9,'HalfPowerFrequency2', ...
                       50.1,'DesignMethod','butter','SampleRate',fs_p);  
    pres(2,nn) = {filtfilt(net_filter,pres{1,nn})} ; 

    % Calibration
    pres(2,nn) = {((1.13*(pres{2,nn}))-31.4)/au} ;

    % Lowpass filter
    [b,a]=butter(4,f_filt/(0.5*fs_p));  
    pres(2,nn) = {filtfilt(b,a,pres{2,nn})};

    % Downsampling
    % pres(2,nn) = {downsample(pres{2,nn}, ds_fact)};    

    % Ask what to display
    answer = questdlg(strcat('Do you want to display pressure channel ' ...
        , num2str(nn),' ?'), ...
        'Select pressure channels', ...
        'Yes','No','No');
    switch answer
        case 'Yes'
            disp(strcat('pressure channel ', num2str(nn),' will be displayed.'))
        case 'No'
            pres(:,nn) = [] ;  % if is is the last cell line, they will ve removed 
            disp(strcat('pressure channel ', num2str(nn),' will NOT be displayed.'))
    end 
end 

