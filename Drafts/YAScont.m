% Detect contractions are calculate outcome measures per conctraction. 
% INPUT
%   - pressure: pressure signal
%   - fs: sampling frequency
%   - name: name of analyzed signal
% OUTPUT = int_cont 
clearvars -except pres t fs 
%%
for i = 1:size(pres,2)
    mod = pres{2,i} ; 
    
    % Find all peaks (contractions) in the pressure signal
    p_peaks = smoothdata(mod,'SmoothingFactor', 0.01);
    [pks_max, locs_max] = findpeaks(p_peaks, fs, 'MinPeakDistance', 5, 'MinPeakProminence',1);
    
     % Find the valleys to determine the start and end times 
    pressure_valleys = smoothdata(mod,'SmoothingFactor', 0.0005);
    [pks_min, locs_min] = findpeaks(-pressure_valleys, fs, 'MinPeakDistance', 4, 'MinPeakProminence', 0.3);  
    
    % Remove peaks that occur in first and last second of signal
    a = find(locs_max < 1);
    b = find(locs_max > (t(end) - 1));
    locs_max(a) = [];
    locs_max(b) = [];
    pks_max(a) = [];
    pks_max(b) = [];
    
    % Determine start times per contractions (based on the found peaks).
    startcont = zeros(length(pks_max),1);
    for ii=1:length(pks_max)
       if sum(locs_min < locs_max(ii)) ~= 0 
           low_locs = locs_min(locs_min < locs_max(ii));
           [~,idx]=min(abs(low_locs-locs_max(ii)));
           minVal=low_locs(idx);
           startcont(ii) = minVal ;
       else 
           locs_max(ii) = NaN ; 
       end  
       endcont = rmmissing(locs_max) ; 
       startcont = nonzeros(startcont) ; 
    end
    int_cont = {[startcont endcont]} ; 
end 
