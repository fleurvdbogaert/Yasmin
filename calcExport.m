function [T,Ts,Tav_pres,Tav_stim] = calcExport(stim,pres,name)

%% Get outcomes 
for i = 1:size(stim,2) % amount of channels used stimulation 
    smod = stim{2,i} ; 
    sfs = stim{3,i} ; 
    time = stim{4,i} ;
    smat = stim{5,i} ; 
    %% Create plot 
    figure(i)
    hold on; 

    t = linspace(1,size(smod,1),size(smod,1)); 
if ~isempty(smat)
    dir = strcat(pwd,'/Calculations') ; 
    file = name{1,1}; 
    filename = sprintf('%s/Stim%d%s%s.xlsx', dir,i,file,datestr(now,'yyyymmdd_HH_MM_ss')); % state filename of the output file
         
    for ii = 1:size(pres,2) % amount of channels used pressure  
        pmod = pres{2,ii} ; 
        pfs = pres{3,ii} ;
        tme = pres{4,ii} ;
        pmat = pres{5,ii} ; 

        tp = (0:numel(pmod)-1)/(pfs); % Time vector 

    if ~isempty(pmat)==1
        hold on ; 
        subplot(size(pres,2),1,ii) ; xlim([0 size(smod,1)/sfs]) ; 

        str = '#80B3FF';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

        for tt = 1:size(time) 
            patch([time(tt,:)./sfs time(tt,2)./sfs time(tt,1)./sfs] ...
                ,[-20 -20 20 20],color,'Edgecolor',color);%'Color','#80B3FF') ; 
        end 
        st1 = sprintf('START CONTRACTION CHANNEL %d',i); 
        st2 = sprintf('END CONTRACTION CHANNEL %d',i); 

        hold on 
        xline(tme(:,1)./pfs,'g-',st1) ; 
        xline(tme(:,2)./pfs,'r-',st2) ; 
        plot(tp, pmod, 'b-', 'LineWidth', 3 );%,'Color','#80B3FF'); % mod data
        xlabel('Time [s]', 'FontSize', 10);
        hold on
        set(gcf, 'Position',  [200, 200, 1000, 400])      % make a rectangular figure
        
%% contractions within intervals

        start = zeros(size(smat,1),1) ; % size of pressure intervals of only contractions 
        stop = zeros(size(smat,1),1) ;
        inter = zeros(size(smat,1),2) ; 
        condur = cell(1,4);

        for ff = 1:size(smat,1) % all stimulation intervals 
            for f = 1:size(tme,1) % contraction intervals 
                
                strt = 0; 
                stp = 0; 
        
                if (tme(f,1)>= smat(ff,2) && tme(f,1)<= smat(ff,3)) == 1 %als het binen de interval ligt 
                    strt = 1 ; % one value means the start (left) or end (right) of the contraction lies within a stimulation
                end 
        
                if (tme(f,2)>= smat(ff,2) && tme(f,2)<= smat(ff,3)) == 1 %als het binen de interval ligt 
                    stp = 1 ; % one value means the start (left) or end (right) of the contraction lies within a stimulation
                end 
                    %contract(ff,1) = (round(sum(sum(strt+stp))/2)) ; % automatically rounds up  
        
                    if strt+stp==0 %sum(bet,2)==0  % contraction start and end are not within the ranges  
                        start(f,1) = 1 ; 
                        stop(f,1)= 1 ;
                    elseif strt+stp==1 %sum(bet,2)==1 % contraction start or end is in the interval 
                        start(f,1) = tme(f,1)*strt ;   st = start(f,1); 
                        stop(f,1) = tme(f,2)*stp ; so =stop(f,1) ; 
        
                        st(start(f,1)==0)=smat(f,2) ;
                        so(stop(f,1)==0)=smat(f,3) ; 
        
                        start(f,1) = st; 
                        stop(f,1) = so ; 

                        inter(f,:) = [tme(f,1) tme(f,2)]; 
        
                    elseif strt+stp==2 %sum(between,2)==2 % contraction start and end are within the interval 
                        stop(f,1) = tme(f,2) ; 
                        start(f,1) = tme(f,1) ; 
                        
                        inter(f,:) = [tme(f,1) tme(f,2)];
                    end 
                if ~isempty(tme)
                    condur(1,1) = {[start stop]} ; % start and end times within stim interval 
                    condur(1,2) = {[tme(f,1) tme(f,2)]} ; % start and end times of the whole contraction
                    condur(1,3) = {stop-start} ; % duration of the contraction within interval 
                    condur(1,4) = {inter(:,2)-inter(:,1)} ; %duration of whole contraction 
                else 
                    condur = [] ; 
                end
            end 
        end 
        %% Individual contractions - Outcomes 
        Xp = {'no contraction', 'contraction'} ;
        Xs = {'no stimulation', 'stimulation'} ;

        % Labels of intervals. 
        label = Xp(1+pmat(:,1))' ;  % 1. labels
        
        % Start time, end time and duration of the contraction
        starttime = (pmat(:,2)-1)./pfs ; starttime(1,1) = 0;      % 2.start
        peaktime = (pmat(:,3)-1)./pfs ;       % 3.stop
        duration = peaktime-starttime ; % 4.duration
        
        % Start pressure, end pressure and absolute height of
        % contraction, all values are for one contraction. 
        startpres = pmod(round(pmat(:,2),0),1) ; 
        peakpres = pmod(round(pmat(:,3),0),1) ; 
        absheight = peakpres - startpres ; % 5.height
        %relheight = startpres./peakpres ;
        %avheight = (startpres+peakpres)/2 ; 
        
        % Slope of the contraction per contraction
        slope = absheight./duration ;    % 6.slope 
        
        % Integral calculated per contraction
        mask = pmat(:,2) <= t & t <= pmat(:,3) ; %create 7x7200000 array with masks 
        modnew = repmat(abs(pmod),1,size(startpres,1))' ; %create 7x7200000 matrix with mod 
        integral = trapz(modnew.*mask,2) ;   % 7. integral: integrate over rows 
        
        %% Averages of all contractions - outcomes 
        
        numcon = sum(pmat(:,1)) ; % total number of contractions
        numstim = sum(smat(:,1)) ; %total number of stimulations 
        
        avdurtot = (sum(duration.*pmat(:,1)))/numcon ; % average contraction time 
        avstartpres = (sum(startpres.*pmat(:,1)))/numcon ; % average start pressure contraction
        avpeakpres = (sum(peakpres.*pmat(:,1)))/numcon ; 
        avabsheight = (sum(absheight.*pmat(:,1)))/numcon ;
        %avrelheight = (sum(relheight.*pmat(:,1)))/numcon ;
        
        avslope = (sum(slope.*pmat(:,1)))/numcon ;    % 6.slope 
        
        %% All stimulations - outcomes  

        % Labels of intervals. 
        labels = Xs(1+smat(:,1))' ;  % 1. labels
        
        poi = condur{:,1} ; 
        % Start time, end time and duration of the contraction
        stim_starttime = ((poi(:,1)-1)./sfs)' ; starttime(1,1) = 0;      % 2.start
        stim_peaktime = ((poi(:,2)-1)./sfs)' ;       % 3.stop
        stim_duration = (stim_peaktime-stim_starttime)' ; % 4.duration
        
        % Start pressure, end pressure and absolute height of
        % contraction, all values are for one contraction.
        
        stim_startpres = smod(round(poi(:,1),0),1)' ; 
        stim_peakpres = smod(round(poi(:,2),0),1)' ; 
        stim_absheight = (stim_peakpres - stim_startpres)' ; % 5.height
        %stim_relheight = stim_startpres./stim_peakpres ;
        %stim_avheight = (stim_startpres+stim_peakpres)/2 ; 
        
        % Slope of the contraction per contraction
        stim_slope = (stim_absheight./stim_duration)' ;    % 6.slope 
        
        % Integral calculated per contraction
        smask = stim_startpres <= t & t <= stim_peakpres ; %create 7x7200000 array with masks 
        smodnew = repmat(abs(pmod),1,size(stim_startpres,1))' ; %create 7x7200000 matrix with mod 
        stim_integral = (trapz(smodnew.*smask,2))' ;   % 7. integral: integrate over rows 
%% Averages over all stimulations - outcomes  

        tot = 0 ; 
        sta = 0 ; 
        sto = 0 ; 

        poj = condur{:,4} ; 
        for h = 1:size(poi,1) % number of intervals within stimulation
            if (smat(h,1)) == 1 %dit is een stimulatie 
       
                tot = tot + poj(h,1) ; %total duration within stim 
                sta = sta + pmod(round(smat(h,2)),1) ; %total start pres within stim 
                sto = sto + pmod(round(smat(h,3)),1) ; %total end pres within stim           
            end 
        end 
        
        Averages = [{'contractions'};{'contractions during stimulations'}]; 

        avdur_c = (tot / numcon)/pfs ; % average duration of whole con per contraction during stimulation
        avdur_s = (tot / numstim)/pfs ; % average duration of contraction per stimulation
        Duration_average = [avdur_c; avdur_s] ; 

        avstapr_c = sta / numcon ; 
        avstapr_s = sta / numstim ; 
        Startpres_average = [avstapr_c ; avstapr_s]; 
        
        avstopr_c = sto / numcon ; 
        avstopr_s = sto / numstim ; 
        Peakpres_average = [avstopr_c ; avstopr_s]; 
        
        avabsheight_c = (sto-sta)/numcon ; 
        avabsheight_s = (sto-sta)/numstim ; 
        Absheight_average = [avabsheight_c ; avabsheight_s];

        avslope_c = avabsheight_c./avdur_c ; 
        avslope_s = avabsheight_s./avdur_s ; 
        Slope_average = [avslope_c;avslope_s];

       % detect miction
       miction = zeros(size(pmat,1),1);
       for m=1:size(pmat,1)
           if pmat(m,1)==1
                mict = pmod(pmat(m,3):(pmat(m,3)+pfs)); 
                Y = fft(mict);
                L = length(mict);
                P2 = abs(Y/L); % fs = length L 
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                fnew = pfs*(0:(L/2))/L;
                [~, loc] = findpeaks(P1, fnew, 'MinPeakHeight',0.5);
                loc = round(loc);
                freqs = 8:1:17;
                c = sum(ismember(loc,freqs));
               if c > 0
                   miction(m,1) = 1;
               end 
           end 
       end
        %% Tables 
            T =  table(label, ...
                starttime, peaktime, duration,...
                startpres, peakpres, absheight,... %relheight,avheight, ...
                slope, integral, miction) ;
            Ts =  table(labels, ...
                stim_starttime, stim_peaktime, stim_duration,...
                stim_startpres, stim_peakpres, stim_absheight,... %stim_relheight,stim_avheight, ...
                stim_slope, stim_integral) ;
            Tav_pres = table(numcon, numstim, ...
                avdurtot, avstartpres, avpeakpres, avabsheight,... %avrelheight, ...
                avslope); 
            Tav_stim = table(Averages, Duration_average, ...
                Startpres_average, ...
                Peakpres_average, ...
                Absheight_average, ...
                Slope_average) ; 

        rng1 = sprintf('A%d:H%d',size(pmat,1)+3,size(pmat,1)+4) ;
        rng2 = sprintf('A%d:J%d',size(pmat,1)+6,size(pmat,1)+8) ;
        rng3 = sprintf('A%d:AG%d',size(pmat,1)+10,size(pmat,1)+(10+size(poi,1))) ;

        warning('off','MATLAB:xlswrite:AddSheet'); %optional
        writetable(T,filename,'Sheet',ii);
        writetable(Tav_pres,filename,'Sheet',ii,'Range',rng1);
        writetable(Tav_stim,filename,'Sheet',ii,'Range',rng2);
        writetable(Ts,filename,'Sheet',ii,'Range',rng3);

    elseif isempty(pmat)==1
    
    writetable(cell2table({'No pressure measurement'}), ...
        filename,'Sheet',ii,'Range','A1:A2') ;  
    
    end
    end     
end 
end 




