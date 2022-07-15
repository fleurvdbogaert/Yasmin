%% Creeren van gelabelde data set 
% Y. (Yasmin) Ben Azouz 
% 4559843 
% Functional Urology, Internship TM-2 
% July 2022 

%% Getting data and labels 
fs_press = 60000;

% Cell with name and data from second pressure measurement 
ROOT_DIR = uigetdir ; 
d = dir(fullfile(ROOT_DIR)) ; % and return all those matching files
d(1:3) = [] ;                 % remove '...', '...' and '.DS.files'
data = cell(3,numel(d)) ;     % create 3x35 cell array, top names bottom data  
for i=1:numel(d)
    data(1,i) = {d(i).name} ; 
    data(2,i) = {importdata(fullfile(ROOT_DIR,d(i).name)).pressure(:,2)} ; %take the second (anal) pressure 
    data(2,i) = {1.13*(data{2,i})-31.4} ;   % apply correction to values

    % Smooth and find peaks 
    [~, plocs] = findpeaks(smoothdata(data{2,i}, ...
        'SmoothingFactor', 0.01), ...
        fs_press,'MinPeakDistance', 5, 'MinPeakProminence',1);

    % Time vector 
    t_press = (0:numel(data{2,i})-1)/(fs_press);

    % Plot figure and peak points 
    figure; 
    h1 = plot(t_press, data{2,i}, '-b', 'LineWidth', 1);
    set(gcf, 'Position',  [200, 200, 1000, 400])      % Set the size of the figure to make it more of a rectangular 
    hold on; 
    xline(plocs, 'LineWidth',2);                 % Plot the onset of each contraction as a vertical line 
    title('Labelling Pressure Peaks'); 
    ylabel('Pressure [cmH_2O]','FontSize', 10);  
    xlabel('Time [s]', 'FontSize', 10); 
 
    % Change, keep or delete peaks 
    for ii=1:numel(plocs)
        xline(plocs(ii), 'b', 'LineWidth',2)

        answer = questdlg('Is the green line a contraction peak?', ...
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
                title('Select peak of contraction!');
                [xn, ~] = ginput(1); 
                plocs(end+1,:) = xn;
                xline(xn, '-');
                disp('Contraction peak added')
            case 'No'
                disp('No more contraction peaks added')
        end    
    end 
    
    plocs = sort(nonzeros(plocs)) ; % remove zeros and sort ascending 

    % data(3,i) = {plocs} ; % Add peaks to 3rd data cell 

    % Add onsets
    choiceon = 'Yes' ; 
    olocs = zeros(2*numel(plocs),1) ;  % create array twice the size of the peaks for enough space for preallocation
    on = 0 ; 

    while isequal(choiceon,'Yes')
        choiceon = questdlg('Do you want to add an onset?', ...
        'Add contraction onset', ...
        'Yes','No','No');
        % Handle response
        switch choiceon
            case 'Yes' 
                title('Select onset of contraction!');
                [xo, ~] = ginput(1); 
                on = on+1 ; 
                olocs(on,:) = xo;
                xline(xo, '-');
                disp('Contraction onset added')
            case 'No'
                disp('No more contraction onsets added')
        end    
    end 
    olocs = sort(nonzeros(olocs)) ;  % remove zeros and sort ascending 
    % data(4,i) = {olocs} ; % add onset data to cells 

    close 

    % Create ROI
if numel(plocs) + numel(olocs) == 2  
    if plocs(1)>olocs(1)
        if isequal(numel(plocs), numel(olocs)) % here every onset is followed by peak
            roi_lims =  cat(2, olocs, plocs) ; 
        elseif numel(plocs) < numel(olocs)  % here we are missing a peak at the end 
            plocs = [plocs; max(t_press)] ;
            roi_lims =  cat(2, olocs, plocs) ; 
        end 
    elseif olocs(1)>plocs(1)
        if numel(olocs) < numel(plocs)      % here we are missing an onset at begin 
            olocs = [0;olocs] ; 
            roi_lims =  cat(2, olocs, plocs) ;
        elseif isequal(numel(plocs), numel(olocs))  % here peak and onset are missing 
            plocs = [plocs; max(t_press)] ;
            olocs = [0;olocs] ; 
        end 
 elseif numel(plocs) + numel(olocs) < 2 
    if isempty(olocs) == 1 
        olocs = zeros(1,1) ; 
        roi_lims = cat(2,olocs,plocs) ; 
    elseif isempty(plocs) == 1 
        plocs = zeros(1,1) ; 
        roi_lims = cat(2,olocs,plocs) ; 
    end 
    end 
end 

    data(3,i) = {zeros(size(data{2,1}))} ; 
    % create label array,zero is nothing, 1 is roi
    for r = 1:length(roi_lims)
        range = round(roi_lims(r,1)*fs_press+1:roi_lims(r,2)*fs_press+1) ; 
        data{3,i}(range) = 1 ; 
    end
end 
