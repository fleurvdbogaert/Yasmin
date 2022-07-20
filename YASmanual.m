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
%%

for i = 1:size(data,2)   % loop over all measurements included 
    if ~isempty(data{1,i}) % if measurements are included

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
            new = sort(nonzeros(pnt(:,ii))) ; % remove zeros and sort ascending
        end 
    end 
    close 
end 


