function [pres_out] = contDetection(pres_in)
% Y. (Yasmin) Ben Azouz
% July 2022 
% 4559843
%% CONTDETECT: automatically detects contractions  
% INPUT
%   - pressure: pressure signal
% OUTPUT = 
%   - int_cont: interval contractions, start and end times 

pres = pres_in ; 
pres_out = pres_in ; 

for i = 1:size(pres,2) % loop over pressure channels  
    if isempty(pres{1,i})==0 % measurement channel used pressure 
        % Variables 
        mod = pres{2,i} ; 
        fs = pres{3,i} ; 
        t = linspace(1,length(mod),length(mod)) ; 
        
        % Find contraction peaks 
        p_peaks = smoothdata(mod,'SmoothingFactor', 0.01);
        [~, locs_max] = findpeaks(p_peaks(fs:t(end)-fs,1), 'MinPeakDistance', ...       % (1:t(end)-1) to remove peaks in first and last seconds
            5, 'MinPeakProminence',1);
        
        % Find the valleys 
        p_valleys = smoothdata(mod,'SmoothingFactor', 0.0005);
        [~, locs_min] = findpeaks(-p_valleys(fs:t(end)-fs,1), 'MinPeakDistance', ...
            4, 'MinPeakProminence', 0.3);  
        
        % Determine start times per contractions 
        startcont = zeros(length(locs_max),1); %preallocate 
        for ii=1:length(locs_max) 
           if sum(locs_min < locs_max(ii)) ~= 0 
               low_locs = locs_min(locs_min < locs_max(ii));
               [~,idx]=min(abs(low_locs-locs_max(ii)));
               minVal=low_locs(idx);
               startcont(ii) = minVal ;
           else 
               locs_max(ii) = NaN ;
               startcont(ii) = NaN ; 
           end  
           endcont = rmmissing(locs_max) ; 
           startcont = nonzeros(startcont) ; 
        end
        int = {[startcont endcont]} ; 
    elseif isempty(pres{1,i})==1 
        int = {[]} ; 
    else 
        int = {NaN(1,2)};
    end
    pres_out(4,i) = int ;
end 
end

