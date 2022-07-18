%% General Information
% Y. (Yasmin) Ben Azouz 
% 4559843
% July 2022 
% If the code errors/does not run, all calculations can be done manually. 
% All detected start and end points will be plotted for accordance and can be adjusted  

%% Plotting and according stimulations 
t_stim = (0:numel(stim)-1)/(fs_press); % Time vector 

%% Plotting and according Pressures 
for i = 1:size(pres,2)
    if ~isempty(pres{1,n})
        t_press = (0:numel(pres{2,n})-1)/(fs_press); % Time vector 

        % figure 
        figure 
        h1 = plot(t_press, pres{2,n}, '-b', 'LineWidth', 1);
        set(gcf, 'Position',  [200, 200, 1000, 400])      % Set the size of the figure to make it more of a rectangular 
        hold on; 
        %% ONSET en piek plotten 
        xline(plocs, 'LineWidth',2);                 % Plot the onset of each contraction as a vertical line 
        title('Labelling Pressure Peaks'); 
        ylabel('Pressure [cmH_2O]','FontSize', 10);  
        xlabel('Time [s]', 'FontSize', 10); 
    end 
end 

 
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



