%% General Information
% Y. (Yasmin) Ben Azouz 
% 4559843
% July 2022 
%% Figures 
% %% Variables
% FS_POT = 60000;        % [Hz] sampling frequency
% FS_PRESS = 60000;
% 
% %% plots 
% figure(1)
% plot(xa,a,...
%     xa(TF),a(TF),'rx') ; 
% 
% figure(2)
% k = 28; 
% x = linspace(1,length(data{2,k}),length(data{2,k})) ; 
% xz = linspace(1,length(data{3,k}),length(data{3,k})) ; 
% figure
% subplot(2,1,1) ; plot(x,data{2,k}) ; title('raw signal')
% subplot(2,1,2) ; plot(xz,data{3,k}) ; title ('movstd')
%% Get stimulation start and end 

ROOT_DIR = uigetdir ; 
d = dir(fullfile(ROOT_DIR)) ; % and return all those matching files
d(1:3) = [] ;                 % remove '...', '...' and '.DS.files'
data = cell(5,numel(d)) ;     % create 4x35 cell array,names,raw,data,peak,start,stim 
for k=1:numel(d)
    data(1,k) = {d(k).name} ; 
    data(2,k) = {importdata(fullfile(ROOT_DIR,d(k).name)).stimulation(:,1)} ; %take first stimulation

    absol = abs(data{2,k}) ; 

    m = movstd(absol,10000) ;
    data(3,k) = {m} ; 
    data(4,k) = {std(m)} ; 

     [~,locs] = findpeaks(data{3,k}, ...
        'MinPeakProminence', (mean(data{3,k}))*3);

    % Remove zeros, they mess with the ranges 
    if isempty(locs) == 0 
    [zks,~] = find(~data{3,k}(1:locs(1))) ;
    n = length(zks) ; 
    div = (n*(n+1))/2 ; 
    if zks ~= 0 
        if isequal(sum(zks), div) 
            data(3,k) = {data{3,k}(locs(1):end)} ;
        end
    end 
        % end of measurement not correct 
    [fks,~] = find(~data{3,k}(locs(end):end)) ;
    f = length(data{3,k}) ; 
    g = length(data{3,k})-length(fks); 
    div = (f*(f+1))/2 - (g*(g+1))/2 ; 
    if fks ~= 0 
        if isequal(sum(fks), div) 
            data(3,k) = {data{3,k}(1:locs(end))} ;
        end 
    end 
        % re-calculate peaks >> should not find peak as half is removed 
    [~,locs] = findpeaks(data{3,k}, ...
        'MinPeakProminence',(mean(data{3,k})*3));
        % re-calculate range 
    data(4,k) = {std(m)} ; %{max(data{3,k})-min(data{3,k})} ;
    end 
%end 

    % STIM OR NOSTIM
    if data{4,k} < 10 % boundary for nostim range 

        if isempty(locs) == 1 
            startstim = [] ; 
            endstim = [] ; 
        elseif isempty(locs) == 0 % waarom zou no stim nog een piek hebben? 
            startstim = [] ; 
            endstim = [] ;
        end 

    elseif data{4,k} > 10 % stim present 

        if isempty(locs) == 1 % cont stim 
            startstim = 1 ;    % beginning stim is one 
            endstim = length(data{3,k}) ;  % end stim is end 

        elseif isempty(locs) == 0
            if (rem(length(locs),2) ~= 0)==1
                % een start heeft links lage std rechts hoge std 
                % een eind heeft links hoge std links lage std 
                startstim = zeros((length(locs)+1)/2) ; 
                endstim = zeros((length(locs)+1)/2) ;

                for l = 1:length(locs)
                    % Take range of 1000? 
                    left = abs(data{2,k}(locs(l)-1000:locs(l))) ; 
                    right = abs(data{2,k}(locs(l):locs(l)+1000)) ; 
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
                    endstim(:,end) = length(data{3,k}) ; 
                elseif sum(startstim==0)~=0 
                    sort(startstim) ; 
                end 
        % even locs > stim started and ended 
            elseif (rem(length(locs),2) == 0) ==1
                startstim = zeros(length(locs)/2) ; 
                endstim = zeros(length(locs)/2) ;

                for l = 1:length(locs)
                    % Take range of 500? 
                    left = abs(data{2,k}(locs(l)-1000:locs(l))) ; 
                    right = abs(data{2,k}(locs(l):locs(l)+1000)) ; 
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
        data(5,k) = {[startstim; endstim]} ; 
    end 
end 
