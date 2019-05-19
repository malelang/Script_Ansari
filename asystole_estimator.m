%% script

if(~isempty(ppg))
    alarmResult = AsysEstimate(ppg,Fs,N0_d,N_d,recordName,ann_ppg,n_ppg_beats);
elseif(~isempty(abp))
    alarmResult = AsysEstimate(abp,Fs,N0_d,N_d,recordName,ann_abp,n_abp_beats);
else
%     fprintf('No PPG nor ABP, i.e. No Belles code\n');
    alarm_result = 1;
end


