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
% mkdir 'Figures' % niet nodig, deze figuren worden niet bekeken 
mkdir 'Calculations'
FS_S = 60000;        % [Hz] sampling frequency stimulation
FS_P = 60000;      % [Hz] sampling frequency pressure  

% pres(2,nn) = {((1.13*(pres{1,nn}))-31.4)/au} ;

% DS_FACT = 2;            % Downsampling factor
F_FILT = 18;            % [Hz] low-pass filter frequency for pressure 

%% Execute 
[STIM,PRES,DIR_OUT] = loadModify(F_FILT,FS_P,FS_S) ;
%%
[INT_STIM] = stimDetection(STIM) ; 
%% 
[INT_PRES] = contDetection(PRES) ; 
%% 
[CHCK_STIM] = manualCheck(INT_STIM,'stimulation') ; 
%%
[CHCK_PRES] = manualCheck(INT_PRES,'pressure') ; 
%% 
[OUT] = calcExport(CHCK_PRES, 'pressure') ; 