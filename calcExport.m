function [data_out] = calcExport(data_in, tpe)
 
data = data_in ; 
data_out = data ; 

for i = 1:size(data,2) % amount of channels used 
    raw = data{1,i} ; 
    mod = data{2,i} ; 
    fs = data{3,i} ; 
    time = data{4,i} ; 

    t = linspace(1,size(mod,1),size(mod,1)); 

    if ~isempty(raw) % measurement channel used

        %% Get intervals 
        num = size(time,1); 
        % Intervals
        if num > 1
            if time(1,1) == 0                   % start activation
                startint = time(:,2)+1 ;        % act end = int start 
                endint = vertcat(time(2:end,1)-1,t(end)) ; % act start = int end  
                if time(end) == t(end)          % end with act 
                    startint = time(1:end-1,2)+1 ; % act end = int start 
                    endint = time(2:end,1)-1 ;  % act start = int end 
                end 
            elseif time(end) == t(end)          % end with act
                startint = vertcat(1, time(1:end-1,2)+1) ; % act end = int start 
                endint = time(:,1)-1 ;          % act start = int end  
            else % no partial activation
                startint = vertcat(1, time(:,2)+1) ; % act end = int start 
                endint = vertcat(time(:,1)-1, t(end)) ; % act start = int end              
            end 
        elseif num == 1 % single %(:,x) should not be necessary really >> add error? In case it is larger?
            if time(1,1) == 1 % partial at start 
                if (time(end)==t(end))==1 || sum(sum(time))==1 % full or none
                   startint = time(:,1); 
                   endint = time(:,2) ; 
                else 
                   startint = time(:,2)+1 ; 
                   endint = t(end) ; 
                end 
            elseif time(1,2) == t(end) % partial at end 
                startint = 1 ; 
                endint = time(:,1)-1 ; 
            else % no partial
                startint = vertcat(1, time(:,2)-1) ; 
                endint = vertcat(time(:,1)+1,t(end)) ; 
            end 
        end 
        %% Outcomes 
        act = horzcat(ones(size(time,1),1),time) ;
        int = horzcat(zeros(size(startint,1),1), startint,endint) ;
        
        ca = vertcat(act,int) ; 
        [~,idx] = sort(ca(:,2));            % sort just the first column
        mat = ca(idx,:);                    % sort the whole matrix using the sort indices
        
        start = mat(:,2)/fs ; start(1) = 0;      % 2.start
        stop = mat(:,3)/fs ;       % 3.stop
        duration = stop-start/fs ; % 4.duration
        height = mod(round(mat(:,3)),1) - mod(round(mat(:,2)),1) ; % 5.height  
        slope = height./duration ;    % 6.slope 
        
        mask = start <= t & t <= stop ; %create 7x7200000 array with masks 
        modnew = repmat(mod,1,size(start,1))' ; %create 7x7200000 matrix with mod 
        integral = trapz(modnew.*mask,2) ;   % 7. integral: integrate over rows 
    
        if isequal(tpe,'pressure')
        
            X = {'interval', 'stimulation'} ; %insert label names press 
            label = X(1+mat(:,1))' ;    % 1. labels pressure 
        
            T =  table(label,start,stop, duration, height, slope, integral) ;  
        elseif isequal(tpe,'stimulation')
            X = {'interval', 'contraction'} ; 
            label = X(1+mat(:,1))' ;  % 1. labels stimulation
            
            % detect miction
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
        
            T =  table(label,start,stop, duration, height, slope, integral, miction) ; 
        else 
            % error message? you have not filled out the right tpe 
        end 

    elseif isempty(raw)
        T = [] ; 
    else 
        %error message? your raw data i contains an abnormal value. Klick
        %... in the workspace, then click the first row of the ith column.
        %What do you see? 
    end 
    data_out(5,i) = {T} ; 
end 



