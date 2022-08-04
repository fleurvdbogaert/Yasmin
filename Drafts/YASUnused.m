%% 04.08.2022 

%                 if f <= size(tme,1) % located on the left side which means start
%                     between(f,1) = tme(f)/fs ; start(1) = 0; % eind stim - start con
%                     between(f,2) = mat(ff,3)/fs ; 
% 
%                 elseif f > size(tme,1) % located on the right side which means contraction end  
%                     between(f,2) = tme(f)/fs ; 
%                     between(f,1) =  mat(ff,2)/fs ; start(1) = 0 ;  % eind con - start stim 
%                 end 
%                 condur(ff,1) = stop-start ; 
%                 
%                 mask = start <= t & t <= stop ; %create 7x7200000 array with masks 
%                 modnew = repmat(abs(mod),1,size(start,1))' ; %create 7x7200000 matrix with mod 
%                 durint = trapz(modnew.*mask,2) ;   % 7. integral: integrate over rows 

%% Plotting and according stimulations LOS - date: 18.07.2022 
t_stim = (0:numel(stim)-1)/(fs_stim); % Time vector >> put inside cell with data. 

% slocs = start points 
% elocs = end points 

%% Plotting and according Pressures LOS - date: 18.07.2022 

% plocs = peak points 
% olocs = onset points 

for i = 1:size(pres,2) % loop over all pressure measurements included 
    if ~isempty(pres{1,n}) % if the pressure measurements are included
        t_press = (0:numel(pres{2,n})-1)/(fs_press); % Time vector 

        % figure 
        figure 
        p1 = plot(t_press, pres{2,n}, '-b', 'LineWidth', 1);
        set(gcf, 'Position',  [200, 200, 1000, 400])      % make a rectangular figure  
        title('Labelling Pressure Peaks'); 
        ylabel('Pressure [cmH_2O]','FontSize', 10);  
        xlabel('Time [s]', 'FontSize', 10);

        hold(p1,'on'); % hold state on 
        %% pieken plotten 
        x1 = xline(plocs, 'LineWidth',2);   
        hold(x1,'on'); % hold state on 

        % Change, keep or delete peaks 
        for ii=1:numel(plocs)
            xline(plocs(ii), 'b', 'LineWidth',2)
    
            answer = questdlg('Is the blue line a contraction peak?', ...
                'Check contractions', ...
                'Keep','Change', 'Delete', 'Keep');
            % Handle response
            switch answer
                case 'Keep'
                    xline(plocs(ii), 'g', 'LineWidth', 2)
                case 'Change'
                    xline(plocs(ii), 'r', 'LineWidth',2)
                    [plocs(ii), ~] = ginput(1); 
                    xline(plocs(ii), 'g', 'LineWidth', 2)
                    disp('Changed contraction')
                case 'Delete'
                    xline(plocs(ii), 'r', 'LineWidth', 2)
                    plocs(ii) = 0;
            end    
        end 

        % Add extra peaks 
        choicep = 'Yes' ; 
        
        while isequal(choicep,'Yes')
            choicep = questdlg('Do you want to add another peak?', ...
            'Add contraction peak', ...
            'Yes','No','No');
            % Handle response
            switch choicep
                case 'Yes'
                    title('Select peak of contraction');
                    [xn, ~] = ginput(1); 
                    plocs(end+1,:) = xn;
                    xline(xn, '-');
                    disp('Contraction peak added')
                case 'No'
                    disp('No more contraction peaks added')
            end    
        end 
    
        plocs = sort(nonzeros(plocs)) ; % remove zeros and sort ascending 
        hold(x1,'off') ; % hold state peak lines off 
        %% onsets plotten
        x2 = xline(olocs, 'LineWidth',2); % detected onset lines 
        x3 = xline(plocs, 'g--', 'LineWidth',2) ; % new peak lines 
        hold([x2 x3],'on'); % hold state on 

        % Change, keep or delete peaks 
        for ii=1:numel(plocs)
            xline(plocs(ii), 'b', 'LineWidth',2)
    
            answer = questdlg('Is the blue line a contraction onset?', ...
                'Check contractions', ...
                'Keep','Change', 'Delete', 'Keep');
            % Handle response
            switch answer
                case 'Keep'
                    xline(plocs(ii), 'g', 'LineWidth', 2)
                case 'Change'
                    xline(plocs(ii), 'r', 'LineWidth',2)
                    [plocs(ii), ~] = ginput(1); 
                    xline(plocs(ii), 'g', 'LineWidth', 2)
                    disp('Changed contraction')
                case 'Delete'
                    xline(plocs(ii), 'r', 'LineWidth', 2)
                    plocs(ii) = 0;
            end    
        end

        % Add onsets
        choiceon = 'Yes' ; 
    
        while isequal(choiceon,'Yes')
            choiceon = questdlg('Do you want to add another onset?', ...
            'Add contraction onset', ...
            'Yes','No','No');
            % Handle response
            switch choiceon
                case 'Yes' 
                    title('Select onset of contraction');
                    [xo, ~] = ginput(1); 
                    olocs(end+1,:) = xo;
                    xline(xo, '-');
                    disp('Contraction onset added')
                case 'No'
                    disp('No more contraction onsets added')
            end    
        end 
        olocs = sort(nonzeros(olocs)) ;  % remove zeros and sort ascending 
    end 
end 

%% FIRST TRY STRUCTURE - date: 17.07.2022 
% findpeaks of moving std 
[~,locs] = findpeaks(data{3,k},'MinPeakProminence',30); % 30 is (denk ik?) goede threshold voor NoStim
 
%% Options A: No peaks found, constant standard deviation so no changes 
% determine start and end stim
if isempty(locs) == 1 
    % differentiate between stim and nostim through range 
    if data{4,k} < 10 % boundary for nostim range 
        startstim = [] ; 
        endstim = [] ; 
        % no stimulation given determined
        % data(5,k) = data(3,k) ;  
    elseif data{4,k} > 10 % boundary for full stim? 
        startstim = 1 ;    % beginning stim is one 
        endstim = length(data{3,k}) ;  % end stim is end 
    end 

%% Option B: peaks found 
elseif isempty(locs) == 0 
    %% ONE: remove not measured (usually one peak but can also be more, check begin and end) 

    % make sure to base of of ONE simulation 
    % removing section where no measurement was done
    % romove this section in all press and stim measurements 
    % start of measurement not correct 
    [zks,~] = find(~data{3,k}(1:locs(1))) ;
        if issorted(zks) == 1 
            data(3,k) = {data{3,k}(locs(1):end)} ;
        end
    % end of measurement not correct 
    [fks,~] = find(~data{3,k}(locs(end):end)) ;
        if issorted(fks,'monotonic') == 1 
            data(3,k) = {data{3,k}(1:locs(end))} ;
        end 
    % re-calculate peaks >> should not find peak as half is removed 
    [~,locs] = findpeaks(data{3,k},'MinPeakProminence',30);
    % re-calculate range 
    data(4,k) = {max(data{3,k})-min(data{3,k})} ;
    
    % Data can now have no more peaks and be no stim or constant stim

    %% TWO: determine start and end stim 
    % odd locs > one stim was not finished 
        % stimulated in the direction of the higher std 
        % take one or end as start or end 
        % maar kan ook 1 zijn met nog een stim er in 
        if rem(locs,2) ~= 0 
            % een start heeft links lage std rechts hoge std 
            % een eind heeft links hoge std links lage std 
            for l = 1:length(locs)
                % Take range of 5000? 
                left = data{3,k}(locs(l)-5000:locs(l)) ; 
                right = data{3,k}(locs(l):locs(l)+5000) ; 
                for ll = 1:length(locs)/2
                    startstim = zeros((length(locs)+1)/2) ; 
                    endstim = zeros((length(locs)+1)/2) ;
                    if left < right %start 
                        startstim(ll) = locs(l) ; 
                    elseif right < left %end 
                        endstim(ll) = locs(l) ; 
                    end 
                end 
            end 

            % find if start or end had zero >> cannot be both (odd)
            % missing start is order in ascending order, zero front 
            % missing end is replace zero with last indice 

            if sum(endstim==0)~=0 
                endstim(:,end) = length(data{3,k}) ; 
            elseif sum(startstim==0)~=0 
                sort(startstim)
            end 
        % even locs > stim started and ended 
        elseif rem(locs,2) == 0 
            for l = 1:length(locs)
                % Take range of 5000? 
                left = data{3,k}(locs(l)-5000:locs(l)) ; 
                right = data{3,k}(locs(l):locs(l)+5000) ; 
                for ll = 1:length(locs)/2
                    startstim = zeros(length(locs)/2) ; 
                    endstim = zeros(length(locs)/2) ;
                    if left < right %start 
                        startstim(ll) = locs(l) ; 
                    elseif right < left %end 
                        endstim(ll) = locs(l) ; 
                    end 
                end 
            end 

        % if locs are again empty now
        elseif isempty(locs) == 1  
        % differentiate between stim and nostim through range 
            if data{4,k} < 10 % boundary for nostim range 
                startstim = [] ; 
                endstim = [] ; 
                % no stimulation given determined
                % What does this mean for the rest of the calculations? 
                % Do we leave the stim start and stop empty? 
            elseif data{4,k} > 10 % boundary for full stim? 
                startstim = 1 ;    % beginning stim is one 
                endstim = length(data{3,k}) ;  % end stim is end
            end 
        end 

end 
data(5,k) = {[startstim; endstim]} ; 
%end  
% x = linspace(1,length(data{2,k}),length(data{2,k})) ; 
%     figure 
% plot(x,data{3,k},...
%     x(locs),data{3,k}(locs),'rx') ; 
%%
%     [zks,~] = find(~data{3,k}(1:locs(1))) ;
%     [fks,~] = find(~data{3,k}(locs(end):end)) ;

%% Garbage from Yasstim - date: 15.07.2022
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

%% Creating labelled set - date: 12.07.2022 
% niet nodig, een normale label array is meer dan genoeg 

% dContraction = signalLabelDefinition('Contraction', ...
%     LabelType='roi') ;
% 
% dStart = signalLabelDefinition('Start', ...
%     LabelType='point') ; 
% dPeak = signalLabelDefinition('Peak', ...
%     LabelType='point') ; 
% 
% dContraction.Sublabels(1) = dStart ; 
% dContraction.Sublabels(2) = dPeak ; 
% 
% 
% lbldefs = dContraction;
% %%
% lss = labeledSignalSet(data(2,1),lbldefs) ; 
%% TVD algorithm - date: 28.04.2022
[Smod, cost] = tvd_mm(potential1_filtered, 2, 4) ; 

figure(2)
x2 = linspace(1, 4, 4) ; 
subplot(2,1,1) ; plot(x2,cost) ; title('cost')
subplot(2,1,2) ; plot(x,Smod) ; title ('stimulation')

%% getTimes lengte ipt proberen aan te passen - date: 20.06.2022
        len = length(2*blocks) - length(ipt) ; 
        add = rand(len,1) ; 
        ipt_add = [ipt add] ; 
        for i=1:2:2*blocks
            tstart  = ts(ipt_add(i));
            tend    = ts(ipt_add(i+1));
            stim = potential(fs_pot*tstart:fs_pot*tend);  
            % Use a fft in calcFreq to find the frequency in the stimulation
            % interval.
            freq = calcFreqV1(stim, fs_pot);
            stim_block(i,:) = [tstart tend];
        end 
