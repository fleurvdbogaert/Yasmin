%% General Information
% Y. (Yasmin) Ben Azouz 
% 4559843
% July 2022 
% >>>> If the code errors/does not run, all calculations can be done manually. 
% All detected start and end points will be plotted for accordance and can be adjusted

%% Plotting and according: general 
%data = cell(4,2) ; % shape of input data 
tpe = 'pressure' ; 
% tpe = pressure or stimulation >> add of axis and titles? 
data = nandata ; 
%%

for i = 1:size(data,2)   % loop over all measurements included 

        raw = data{1,i} ; % raw data 
        mod = data{2,i} ; % modified data 
        pnt = data{3,i} ; % points
        fs = data{4,1} ; % put sample frequency inside cell with data.

        t = (0:numel(raw)-1)/(fs); % Time vector

        % figure  
           % plot figure 
            figure 
            plot(t, mod, '-b', 'LineWidth', 1); % mod data 
            hold on
            set(gcf, 'Position',  [200, 200, 1000, 400])      % make a rectangular figure
    if sum(isnan(pnt)) > 1  % if the data did not go through the stim or cont detection
        for ii = 1:size(pnt,2) % hoort eigenlijk altijd 2 te zijn, gaat om begin en eind. 
    
               % determine name for pop-ups 
               if (isequal('stimulation',tpe)) == 1 
                    ylabel('Stimulation Voltage [AU]','FontSize', 10); 
                    if ii == 1 
                        name = ' STIMULATION START' ; 
                    elseif ii == 2 
                        name = ' STIMULATION END' ; 
                    end 
                elseif (isequal('pressure',tpe)) == 1
                    ylabel('Pressure [cmH_2O]','FontSize', 10);
        
                    if ii == 1 
                        name = ' CONTRACTION ONSET' ; 
                    elseif ii == 2 
                        name = ' CONTRACTION PEAK' ; 
                    end 
               end 
          
            % Manual adding  
            choicep = 'Yes' ; 
            
            while isequal(choicep,'Yes')
                choicep = questdlg(strcat('Do you want to add a ' ...
                    ,name,'?'), ...
                    strcat('Add',name), ...
                'Yes','No','No');
                % Handle response
                switch choicep
                    case 'Yes'
                        title(strcat('Select',name));
                        [xn, ~] = ginput(1); 
                        pnt(end+1,ii) = xn;
                        xline(xn, '-');
                        disp(strcat('Added',name)) 
                    case 'No'
                        disp(strcat('No more',name,' added.'))
                end    
            end 
            hold on ; 
            pnt_left = sort(rmmissing(nonzeros(pnt(:,1)))) ; % remove zeros and sort ascending
            pnt_right = sort(rmmissing(nonzeros(pnt(:,2)))) ;
        end 
            
        
    elseif isempty(raw) && sum(sum(isnan(pnt))) == 0 % if measurements are included (all is filled with [] if not included) 
        for ii = 1:size(pnt,2) % hoort eigenlijk altijd 2 te zijn, gaat om begin en eind. 
               % determine name for pop-ups 
               if (isequal('stimulation',tpe)) == 1 
                    ylabel('Stimulation Voltage [AU]','FontSize', 10); 
                    if ii == 1 
                        name = ' STIMULATION START' ; 
                    elseif ii == 2 
                        name = ' STIMULATION END' ; 
                    end 
                elseif (isequal('pressure',tpe)) == 1
                    ylabel('Pressure [cmH_2O]','FontSize', 10);
        
                    if ii == 1 
                        name = ' CONTRACTION ONSET' ; 
                    elseif ii == 2 
                        name = ' CONTRACTION PEAK' ; 
                    end 
               end 
                
                xline((pnt(:,ii)/fs), 'LineWidth',2);   % lines 
                % xline((pnt(:,2)/fs), 'LineWidth',2);   % lines 
                hold on; % hold state on 
    
                title(strcat('Labelling',name)); 
                   
                xlabel('Time [s]', 'FontSize', 10);
    
                % Change, keep or delete  
                for iii=1:numel(pnt(:,ii))
                    xline((pnt(iii,ii)/fs), 'b', 'LineWidth',2)
            
                    answer = questdlg(strcat('Is the blue line a' ...
                        ,name,'?'), ...
                        'Change, keep or delete?', ...
                        'Keep','Change', 'Delete', 'Keep');
                    % Handle response
                    switch answer
                        case 'Keep'
                            xline((pnt(iii,ii)/fs), 'g', 'LineWidth', 2)
                            disp(strcat('Kept',name))
                        case 'Change'
                            xline((pnt(iii,ii)/fs), 'r', 'LineWidth',2)
                            [pnt(iii,ii), ~] = ginput(1); 
                            xline((pnt(iii,ii)/fs), 'g', 'LineWidth', 2)
                            disp(strcat('Changed',name))
                        case 'Delete'
                            xline((pnt(iii,ii)/fs), 'r', 'LineWidth', 2)
                            pnt(iii,ii) = 0;
                            disp(strcat('Deleted',name))
                    end    
                end 
        
                % Add extra  
                choicep = 'Yes' ; 
                
                while isequal(choicep,'Yes')
                    choicep = questdlg(strcat('Do you want to add another' ...
                        ,name,'?'), ...
                        'Add extra', ...
                    'Yes','No','No');
                    % Handle response
                    switch choicep
                        case 'Yes'
                            title(strcat('Select',name));
                            [xn, ~] = ginput(1); 
                            pnt(end+1,ii) = xn;
                            xline(xn, '-');
                            disp(strcat('Added',name)) 
                        case 'No'
                            disp(strcat('No more',name,' added.'))
                    end    
                end 
                hold on ; 
                pnt_left = sort(rmmissing(nonzeros(pnt(:,1)))) ; % remove zeros and sort ascending
                pnt_right = sort(rmmissing(nonzeros(pnt(:,2)))) ;
        end 
    end 
        %else 
    
        if isempty(pnt_left) == 1 
            pnt_left = zeros(1,1) ; 
            man = cat(2,pnt_left,pnt_right) ; 
        elseif isempty(pnt_right) == 1 
            pnt_right = zeros(1,1) ; 
            man = cat(2,pnt_left,pnt_right) ; 
        else

            if pnt_right(1)>pnt_left(1)
                if isequal(numel(pnt_right), numel(pnt_left)) % here every onset is followed by peak
                    man =  cat(2, pnt_left, pnt_right) ; 
                elseif numel(pnt_right) < numel(pnt_left)  % here we are missing a peak at the end 
                    pnt_right = [pnt_right; max(t_press)] ;
                    man =  cat(2, pnt_left, pnt_right) ; 
                end 
    
            elseif pnt_left(1)>pnt_right(1)
                if numel(pnt_left) < numel(pnt_right)      % here we are missing an onset at begin 
                    pnt_left = [0;pnt_left] ; 
                    man =  cat(2, pnt_left, pnt_right) ;
                elseif isequal(numel(pnt_right), numel(pnt_left))  % here peak and onset are missing 
                    pnt_right = [pnt_right; max(t_press)] ;
                    pnt_left = [0;pnt_left] ; 
                end 
            end 

        end             
   
 
    close 
end 
 


