function [T,Tav_pres,Tav_stim] = calcExport(stim,pres)
  
% data_out = data ; 

%% Get outcomes 
for i = 1:size(stim,2) % amount of channels used stimulation
%     raw = data{1,i} ; 
%     mdd = data{2,i} ;  
%     fs = data{3,i} ; 
%     time = data{4,i} ; 
 
    sraw = stim{1,i} ; 
    smod = stim{2,i} ; 
    sfs = stim{3,i} ; 
    time = stim{4,i} ;
    smat = stim{5,i} ; 

    t = linspace(1,size(smod,1),size(smod,1)); 
    tables = cell(1,size(pres,2)) ;
         
    for ii = 1:size(pres,2) % amount of channels used pressure 
        praw = pres{1,ii} ; 
        pmod = pres{2,ii} ; 
        pfs = pres{3,ii} ;
        tme = pres{4,ii} ;
        pmat = pres{5,i} ; 

%         if ~isempty(mod)
%             
%  
%         elseif isempty(mod)
%             T = [] ; 
%             continue 
%         else 
%             continue 
%         end 
%         
%         tables(1,ii) = {T} ; %put a table for each pressure measurement within a cell structure 
% 
%     end 
%     data_out(5,i) = {tables} ; %put cell with tables 
% end 

%% contractions within intervals

start = zeros(size(tme,1),1) ; % size of pressure intervals of only contractions 
stop = zeros(size(tme,1),1) ;

%contract = zeros(size(smat,1),1); % size of ALL the stimulation intervals 
condur = cell(size(smat,1),1);
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
                start(f,1) = 0 ; 
                stop(f,1)= 0 ;
            elseif strt+stp==1 %sum(bet,2)==1 % contraction start or end is in the interval 
                start(f,1) = tme(f,1)*strt ;   st = start(f,1); 
                stop(f,1) = tme(f,2)*stp ; so =stop(f,1) ; 

                st(start==0)=smat(ff,2) ;
                so(stop==0)=smat(ff,3) ; 

                start(f,1) = st; 
                stop(f,1) = so ; 

            elseif strt+stp==2 %sum(between,2)==2 % contraction start and end are within the interval 
                stop(f,1) = tme(f,2) ; 
                start(f,1) = tme(f,1) ; 
                % add values same ad above? 
            end 
        condur(ff,1) = {[start stop]} ; % start and end times within stim interval 
        condur(ff,2) = {[tme(f,1) tme(f,2)]} ; % start and end times of the whole contraction
        condur(ff,3) = {stop-start} ; % duration of the contraction within interval 
        condur(ff,4) = {tme(f,2)-tme(f,1)} ; %duration of whole contraction 
        
        [T,Tav_pres,Tav_stim] = Outcomes(pmod, smat, pmat, condur, pfs) ; 
    end

end 
    end 
end 
end 

%percdur = duration./condur ; 

%filename = sprintf('%s\\Stimulation information %s %s.xlsx', dir, s.name, datestr(now,'yyyy-mm-dd HH-MM-ss')); % state filename of the output file
%filename_p = sprintf('%s\\Pressure calculations %s %s.xlsx', dir,s.name, datestr(now,'yyyy-mm-dd HH-MM-ss'));
%    warning('off','MATLAB:xlswrite:AddSheet'); %optional
%    writetable(M1,'test.xlsx','Sheet',1);
%    writetable(M2,'test.xlsx','Sheet',2);
% end 

%% 



