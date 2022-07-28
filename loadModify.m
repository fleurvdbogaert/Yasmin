function [stim,pres,dir_outp] = loadModify(au, ds_fact, f_filt, fs_press)
% Data is loaded into cells and modified. 
%   Y. (Yasmin) Ben Azouz
%   4559843
%   July 2022 
%   Only one file can be selected at a time 

%% Variables
fs_new = fs_press / ds_fact;

%% Select input file and load data 
[name,path] = uigetfile('*.*') ;
name = {name} ; 
data = importdata(fullfile(path,name{1,1})) ;

% Display selected file in command window
if isequal(name{1,1},0)
   error('File selection canceled') ;
else 
   disp('File selected: '), disp(name{1,1}) ; 
end
%% Select output directory 
dir_outp = uigetdir(dir_export, 'Select OUTPUT directory');

%% Modify data 

% stimulation 
num_stim = size(data.stimulation,2) ; % amount of stimulation measurements
stim = cell(1,num_stim) ;
for n = 1:num_stim
    stim(1,n) = {data.stimulation(:,n)} ; 

    % Ask what to display
    answer = questdlg(strcat('Do you want to display stimulation channel ' ...
        , num2str(n),' ?'), ...
        'Select stimulation channels', ...
        'Yes','No','Yes');
    switch answer
        case 'Yes'
            disp(strcat('stimulation channel ', num2str(n),' will be displayed.'))
        case 'No'
            stim(n) = [] ; 
            disp(strcat('stimulation channel ', num2str(n),' will NOT be displayed.'))
    end   
end 
       
% pressure 
num_pres = size(data.pressure,2) ; % amount of pressure measurements 
pres = cell(2,num_pres) ; % first raw, second modified 
for nn = 1:num_pres
    pres(1,nn) = {data.pressure(:,nn)} ; 

    % Calibration
    pres(2,nn) = {((1.13*(pres{1,nn}))-31.4)/au} ;

    % Lowpass filter
    [b,a]=butter(4,f_filt/(0.5*fs_new));  
    pres(2,nn) = {filtfilt(b,a,pres{2,nn})};

    % Downsampling
    pres(2,nn) = {downsample(pres{2,nn}, ds_fact)};    

    % Ask what to display
    answer = questdlg(strcat('Do you want to display pressure channel ' ...
        , num2str(nn),' ?'), ...
        'Select pressure channels', ...
        'Yes','No','Yes');
    switch answer
        case 'Yes'
            disp(strcat('pressure channel ', num2str(nn),' will be displayed.'))
        case 'No'
            pres(:,nn) = [] ; 
            disp(strcat('pressure channel ', num2str(nn),' will NOT be displayed.'))
    end 
end 

