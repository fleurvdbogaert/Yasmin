


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