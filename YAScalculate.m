%%   Y. (Yasmin) Ben Azouz
%   4559843
%   July 2022 

%% Determine pressures at stimulation interval start and end times 
% Determine the pressure values at start and end times of stimulation
% intervals

for i = 1:size(stim,2) % amount of stimulaion channels used 
    if isempty(stim{1,1})==0 % measurement channel used stimulation 
        sse = stim{3,i} ; %stimulation start and end points 
        stime = round(stim{3,i}*fs) ; %timepoints start and end stimulations back to sample axis 

        OUTnum_stim = size(sse) ; %number of stimulations 

    for ii = 1:size(pres,2) 
        if isempty(pres{1,1})==0 % measurement channel used pressure 

            % Rename variables 
            pse = pres{3,ii} ; % contraction start and end points 
            ptime = round(pse*fs) ; % start and end times contractions to sample axis 
            p_mod = pres{2,ii} ; % modified pressure signal 

            OUTp_ystim = zeros(size(sse)) ; % make an array the same size to fill with heights stim points 
            OUTp_ypres = zeros(size(pse)) ; % make an array the same size to fill with heights pres points
                
                % Sections stimulations 
                if OUTnum_stim > 1 % multiple stimulation blocks 
                    if stime(1,1) == 0 % start with stimulation
                        startint = sort(vertcat( ...
                            stime(:,1),(stime(:,2)+1))) ; 
                        endint = sort(vertcat( ...
                            stime(:,2),(stime(2:end,1)-1),t(end))) ; 

                        if stime(end) == t(end) % end with stimulation
                            startint = sort(vertcat( ...
                                stime(:,1),(stime(1:end-1,2)+1))) ;  
                            endint = sort(vertcat( ...
                                stime(:,2), stime(2:end,1)-1)) ; 
                        end 
                    elseif stime(end) == t(end) % end with stimulation
                        startint = sort(vertcat( ...
                            stime(:,1),(stime(1:end-1,2)+1),0)) ; 
                        endint = sort(vertcat( ...
                            stime(:,2),(stime(:,1)-1))) ; 
                    else % no partial stimulations
                        startint = sort(vertcat( ...
                            0, stime (:,1), (stime(:,2)+1))); % stim end = int start 
                        endint = sort(vertcat( ...
                            stime(:,2),(stime(:,1)-1), t(end))) ; %stim start = int end              
                    end 
                elseif OUTnum_stim == 1 % one stimulation block
                    if stime(1,1) == 0 
                        sect = [stime ; stime(:,2)+1 t(end)] ;
                        if stime(end)==t(end)
                          sect = stime ;
                        end 
                    elseif stime(1,2) == t(end) % partial at end 
                        sect = [0 stime(:,1)-1 ; stime] ; 
                    else % no partial stimulations
                        sect = [0 stime(:,1)-1 ; stime ; stime(:,2)+1 t(end)] ; 
                    end 
                %else? 
                end 
                stim_sect = endint - startint ; % duration

            for s = 1:size(stime,1) % nog gebruiken? de height geeft dit in principe ook aan en daar ga je over heen en het bespaart een loop
                OUTp_ystim(s,1) = p_mod(:,stime(s,1)) ; % start stim value pressure 
                OUTp_ystim(s,2) = p_mod(:,stime(s,2)) ; % end stim value pressure 

            end 
%%%%%%%%%%% 1. Number of contractions %%%%%%%%%%
            OUTnum_cont = size(pse,1) ; % number of contractions detected 
            OUTcont_pm = size(pse,1) / ((size(p_mod,2)/fs)/60); % number of contractions per minute 

%%%%%%%%%%% 2. Duration %%%%%%%%%%

            % Contractions
            OUTp_dur = (ptime(:,2) - ptime(:,1))/fs ; % durations in samplerate divided by fs gives time axis

            % Intervals 
            if OUTnum_cont > 1
                if ptime(1,1) == 0 % start with contraction
                    startint = ptime(:,2) ; % cont end = int start 
                    endint = [ptime(2:end,1) t(end)] ; %cont start = int end  
                    if ptime(end) == t(end) % end with contraction
                        startint = ptime(1:end-1,2) ; % cont end = int start 
                        endint = ptime(2:end,1) ; %cont start = int end 
                    end 
                elseif ptime(end) == t(end) % end with contraction
                    startint = [0 ptime(1:end-1,2)] ; % cont end = int start 
                    endint = ptime(:,1) ; %cont start = int end  
                else % no partial contractions
                    startint = [0 ptime(:,2)] ; % cont end = int start 
                    endint = [ptime(:,1) t(end)] ; %cont start = int end              
                end 
            elseif OUTnum_cont == 1 % single contraction
                if ptime(1,1) == 0 % partial at start 
                    startint = ptime(:,2) ; 
                    endint = t(end) ; 
                    if ptime(1,2) == t(end) 
                       startint = NaN ; 
                       endint = NaN ; 
                    end 
                elseif ptime(1,2) == t(end) % partial at end 
                    startint = 0 ; 
                    endint = ptime(:,1) ; 
                else % no partial contractions 
                    startint = [0 ptime(:,2)] ; 
                    endint = [ptime(:,1) t(end)] ; 
                end 
            %else? 
            end 
            OUTpint_dur = endint - startint ; 

        for k = 1:size(ptime,1) %loop over the start and end times pressure 

%%%%%%%%%%% 3. Height %%%%%%%%%%
            OUTp_ypres(k,1) = p_mod(ptime(k,2)) - p_mod(ptime(k,1)) ; 
%%%%%%%%%%% 4. Slope %%%%%%%%%% 
            OUTp_slope = OUTp_ypres / OUTp_dur ; 
%%%%%%%%%%% 5. Integral %%%%%%%%%% 
            cont = p_mod(ptime(k,1)):p_mod(ptime(k,2)) ; % interval y-as 
            OUTp_tra = trapz(cont) ; 
%%%%%%%%%%% 6. Voiding %%%%%%%%%%
            mict = p_mod(ptime(k,2)):p_mod(ptime(k,2)+fs); 
            Y = fft(mict);
            P2 = abs(Y/fs); % fs = length L 
            P1 = P2(1:fs/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = fs_press*(0:(fs/2))/fs;
            [~, loc] = findpeaks(P1, f, 'MinPeakHeight',0.5);
            loc = round(loc);
            freqs = 8:1:17;
            c = sum(ismember(loc,freqs));
               if c ~= 0
                   miction = 1;
               else 
                   miction = 0;
               end 
        end         
        end 
    end 
    end 
end 
