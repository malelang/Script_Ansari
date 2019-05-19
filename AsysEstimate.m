

function alarmResult = AsysEstimate(ppg,Fs,N0_d,N_d,recordName,ann_ppg,n_ppg_beats)

unusual_HR_flag=0;
alarmResult = 1;
ppg_win = ppg(N0_d:N_d);

peakp = ann_ppg(n_ppg_beats)';

if(length(peakp) > 2)
    peak_dist = diff(peakp);
    peak_dist(end+1) = N_d - ann_ppg(n_ppg_beats(end));
    
    if(max(60*Fs/(peak_dist')) > 140 || max(60*Fs/(peak_dist')) < 20)
%         fprintf('unusual HR %2.2f \n',max(60*Fs/(peak_dist')));
        unusual_HR_flag=1;
    end
    
%     figure; subplot(2,1,1); plot(ppg_win); title(recordName)
%     subplot(2,1,2); plot([diff(peakp),length(ppg_win)-peakp(end)]);
else
%     fprintf('No peaks detected in PPG');
    unusual_HR_flag=1;
    peak_dist = 0;
end

% figure; subplot(2,1,1); plot(ppg_win); title(recordName)
% subplot(2,1,2); plot(peak_dist); 

if(~isempty(find(peak_dist > 4*Fs)) || unusual_HR_flag == 1 ) %finding if there is a gap of 4 sec or more between peaks
%     fprintf('Asystole! PPG\n');
    alarmResult=1;
else
%     fprintf('NOT Asystole! PPG\n');
    alarmResult=0;
end
end
% close all;