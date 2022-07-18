function [int_stim] = stimDetection(stim)
% This functions detects the start and end times of stimulations 

%   stim = double with raw data from a stimulation measurement 
%   int_stim = 2x(amount of stimulations) matrix, start and end times of the detected stimulations 

absol = abs(stim) ;         % absolute value 
mov = movstd(absol,10000) ;   % standard deviation of moving standard deviation
m = std(mov) ; 

[~,locs] = findpeaks(mov, ...
    'MinPeakProminence', (mean(mov))*3);  % detect start and ends 

% Remove zeros, they mess with the ranges 
if isempty(locs) == 0 
[zks,~] = find(~mov(1:locs(1))) ;
n = length(zks) ; 
div = (n*(n+1))/2 ; 
if zks ~= 0 
    if isequal(sum(zks), div) 
        mov = mov(locs(1):end) ;
    end
end 
    % end of measurement not correct 
[fks,~] = find(~mov(locs(end):end)) ;
f = length(mov) ; 
g = length(mov)-length(fks); 
div = (f*(f+1))/2 - (g*(g+1))/2 ; 
if fks ~= 0 
    if isequal(sum(fks), div) 
        mov = mov(1:locs(end)) ;
    end 
end 
    % re-calculate peaks >> should not find peak as half is removed 
[~,locs] = findpeaks(mov, ...
    'MinPeakProminence',(mean(mov)*3));
    % re-calculate range 
m = std(mov) ; 
end 


% STIM OR NOSTIM
if m < 10 % boundary for nostim range 

    if isempty(locs) == 1 
        startstim = [] ; 
        endstim = [] ; 
    elseif isempty(locs) == 0 % waarom zou no stim nog een piek hebben? 
        startstim = [] ; 
        endstim = [] ;
    end 

elseif m > 10 % stim present 

    if isempty(locs) == 1 % cont stim 
        startstim = 1 ;    % beginning stim is one 
        endstim = length(mov) ;  % end stim is end 

    elseif isempty(locs) == 0
        if (rem(length(locs),2) ~= 0)==1
            % een start heeft links lage std rechts hoge std 
            % een eind heeft links hoge std links lage std 
            startstim = zeros((length(locs)+1)/2) ; 
            endstim = zeros((length(locs)+1)/2) ;

            for l = 1:length(locs)
                % Take range of 1000? 
                left = abs(stim(locs(l)-1000:locs(l))) ; 
                right = abs(stim(locs(l):locs(l)+1000)) ; 
                for ll = 1:(length(locs)+1)/2
                    if sum(left) < sum(right) %start 
                        startstim(ll) = locs(l) ; 
                    elseif sum(right) < sum(left) %end 
                        endstim(ll) = locs(l) ; 
                    end 
                end 
            end 
        

        % find if start or end had zero >> cannot be both (odd)
        % missing start is order in ascending order, zero front 
        % missing end is replace zero with last indice 

            if sum(endstim==0)~=0 
                endstim(:,end) = length(mov) ; 
            elseif sum(startstim==0)~=0 
                sort(startstim) ; 
            end 
    % even locs > stim started and ended 
        elseif (rem(length(locs),2) == 0) ==1
            startstim = zeros(length(locs)/2) ; 
            endstim = zeros(length(locs)/2) ;

            for l = 1:length(locs)
                % Take range of 500? 
                left = abs(stim(locs(l)-1000:locs(l))) ; 
                right = abs(stim(locs(l):locs(l)+1000)) ; 
                for ll = 1:length(locs)/2
                    if sum(left) < sum(right) %start 
                        startstim(ll) = locs(l) ; 
                    elseif sum(right) < sum(left) %end 
                        endstim(ll) = locs(l) ; 
                    end 
                end 
            end 
        end 
    end 
    int_stim = {[startstim; endstim]} ; % put inside cell to move 
end 
% end of function
end 


