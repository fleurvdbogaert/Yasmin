function [T,Tav_pres,Tav_stim] = Outcomes(moddata, smatdata, pmatdata, duras, fs)
%OUTCOMES Summary of this function goes here
%   Detailed explanation goes here

% if isequal(tpe,'pressure')
    %X = {'interval', 'contraction'} ; 
% elseif isequal(tpe,'stimulation')
%     X = {'interval', 'stimulation'} ; 
% end 

mod = moddata ; % modified pressure {2,i} array data 
mat = pmatdata ; % mat data {5,i} pressure 
mtt = smatdata ; % mat data {5,i} stimulation
condur = duras ; 

%% Individual contractions - Outcomes 
if ~isempty(mat)==1
% Labels of intervals. 
%label = X(1+mat(:,1))' ;  % 1. labels

% Start time, end time and duration of the contraction
starttime = mat(:,2)./fs ; starttime(1,1) = 0;      % 2.start
stoptime = mat(:,3)./fs ;       % 3.stop
duration = stoptime-starttime ; % 4.duration

% Start pressure, end pressure and absolute height of
% contraction, all values are for one contraction. 
startpres = mod(round(mat(:,2)),1) ; 
stoppres = mod(round(mat(:,3)),1) ; 
absheight = stoppres - startpres ; % 5.height
relheight = startpres./stoppres ;   
avheight = mean(mod(round(mat(:,2):mat(:,3)),1)) ; 

% Slope of the contraction per contraction
slope = absheight./duration ;    % 6.slope 

% Integral calculated per contraction
mask = start <= t & t <= stop ; %create 7x7200000 array with masks 
modnew = repmat(abs(mod),1,size(start,1))' ; %create 7x7200000 matrix with mod 
integral = trapz(modnew.*mask,2) ;   % 7. integral: integrate over rows 

%% Averages of all contractions - outcomes 

numcon = sum(mat(:,1)) ; % total number of contractions
numstim = sum(mtt(:,1)) ; %total number of stimulations 

avdurtot = sum((duration.*mat(:,1)))/numcon ; % average contraction time 
avstartpres = sum(startpres.*mat(:,1))/numcon ; % average start pressure contraction
abstoppres = sum(stoppres.*mat(:,1))/numcon ; 
avabsheight = sum(absheight.*mat(:,1))/numcon ;
avrelheight = sum(relheight.*mat(:,1))/numcon ;

avslope = sum(slope.*mat(:,1))/numcon ;    % 6.slope 

%% Averages over all stimulations - outcomes   
tot = 0 ; 
sta = 0 ; 
sto = 0 ; 

for i = 1:size(condur,1) % number of intervals within stimulation
    if (mtt(i,1)) == 1 %dit is een stimulatie 

        tot = tot + condur{i,4} ; %total duration within stim 
        sta = sta + mod(round(mtt(i,2),1)) ; %total start pres within stim 
        sto = sto + mod(round(mtt(i,3),1)) ; %total end pres within stim
%     elseif (mtt{i,1}) == 0 
%         oto = oto +condur{i,4}
%         avdurint = oto/(size(mtt,1)-num)
    end 
end 

avdur_c= tot / numcon ; % average duration of whole con per contraction during stimulation
avdur_s = tot / numstim ; % average duration of contraction per stimulation

avstapr_c = sta / numcon ; 
avstapr_s = sta / numstim ; 

avstopr_c = sto / numcon ; 
avstopr_s = sto / numstim ; 

avabsheight_c = (sto-sta)/numcon ; 
avabsheight_s = (sto-sta)/numstim ; 

avslope_c = avabsheight_c/avdur_c ; 
avslope_s = avabsheight_s/avdur_s ; 

    T =  table...%(label, ...
        (starttime, stoptime, duration,...
        startpres, stoppres, absheight, relheight, avheight, ...
        slope, integral); %, miction) ;
    Tav_pres = table(numcon, numstim, ...
        avdurtot, avstartpres, abstoppres, avabsheight, avrelheight, ...
        avslope); 
    Tav_stim = table(avdur_c, avdur_s, ...
        avstapr_c, avstapr_s, ...
        avstopr_c, avstopr_s, ...
        avabsheight_c, avabsheight_s, ...
        avslope_c, avslope_s) ; 
elseif isempty(mat)
    T = [] ; 
    Tav_pres = [] ; Tav_stim = [] ; 
end 

%else 
    % error message? you have not filled out the right tpe 
end
%                 % detect miction
%                 mict = p_mod(ptime(k,2)):p_mod(ptime(k,2)+fs); 
%                 Y = fft(mict);
%                 P2 = abs(Y/fs); % fs = length L 
%                 P1 = P2(1:fs/2+1);
%                 P1(2:end-1) = 2*P1(2:end-1);
%                 f = fs_press*(0:(fs/2))/fs;
%                 [~, loc] = findpeaks(P1, f, 'MinPeakHeight',0.5);
%                 loc = round(loc);
%                 freqs = 8:1:17;
%                 c = sum(ismember(loc,freqs));
%                    if c ~= 0
%                        miction = 1;
%                    else 
%                        miction = 0;
%                    end 
