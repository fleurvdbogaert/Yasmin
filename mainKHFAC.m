%% General Information
% Author: Yasmin Ben Azouz, yasminbenazouz@hotmail.com
% Technical Medicine Internship May 2022 - August 2022 

% Previous authors: 
    % Max Ligtenberg, max.ligtenberg@outlook.com, June - Aug 2020
    % Bart Formsma, bartformsma@hotmail.com, Sep - Nov 2020
    % Anne Meester, annemeester95@hotmail.com, Sep - Nov 2021
    % Sabine Josemans, shjosemans@gmail.com, Nov 2021 - Feb 2022
    
% ErasmusMC, dept. Urology, group: functional Urology lab
% Edited in MATLAB R2021b

close all; clear; clc;

%% Variables
% DIR_IMPORT = [];        % when desired, set a default directory to select files from, otherwise use = [];
% DIR_EXPORT = [];        % when desired, set a default directory the OUTPUT directory user interface will open up, otherwise use = [];
mkdir 'Figures'
mkdir 'Calculations'
FS_POT = 60000;        % [Hz] sampling frequency
FS_PRESS = 60000;

%% Input values
C_press1 = 0;         
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
[SData_mod, FS_PRESS] = modifyDataV4(SData, C_press1, C_press2, DS_FACT, FS_PRESS, FS_POT, F_FILT, HF_BLOCKS, LF_BLOCKS, PRESS1, PRESS2);
%%
plotDataV1(SData_mod,FS_POT, FS_PRESS, C_press1, C_press2);
%%
[SData_auto, start_stim, end_stim] = automaticIntervalSelectionV1(SData_mod,FS_POT, FS_PRESS, HF_BLOCKS, LF_BLOCKS, PRESS1, PRESS2);
%%
[SData_exp] = calcOutcomeV4(SData_auto, FS_POT, FS_PRESS, TX, HF_BLOCKS, LF_BLOCKS, PRESS1, PRESS2, start_stim, end_stim);
%%
exportFuncV4(SData_exp, DIR_EXPORT, PRESS1, PRESS2);

msgbox('Operation Completed','Success')
disp('Operation completed!')