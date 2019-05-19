function Phase_Wrapper_batch(alarm)
%% Load Alarms
% If input is provided, only annotate those files; otherwise, annotate
% all.

if(strcmp(alarm,'VF'))
    alarmname = 'Ventricular_Flutter_Fib';
elseif(strcmp(alarm,'VT'))
    alarmname = 'Ventricular_Tachycardia';
end    


alarms = readtable('ALARMS','ReadVariableNames',0,'ReadRowNames',1);    
% if(nargin == 1)
%     alarms = alarms(files,:);
% else
    del = [];
    for i = 1:size(alarms,1)
        if(~strcmp(alarms{i,1},alarmname)) %change the type of alarm here
            del = [del i];
        end
    end
    alarms(del,:) = []; % uncomment if you wanna only annotate the
    % alarms that are specified above

    % alarms = alarms(randperm(size(alarms,1)),:);   % uncomment if you
    % wanna annotate is random order
%end    

% form the list of files, the alarm_types and truth (true or false
% alarm)
files = alarms.Properties.RowNames;
alarm_types = alarms{:,1};
truth = cell(size(files));
for i = 1:length(files)
    if(alarms{i,2}==0)
        truth{i} = 'False';
    else
        truth{i} = 'True';
    end

    alarm_types{i} = strrep(alarm_types{i},'_','\_');
end
%% Run this section for all subjects
phase_wrapped_struct = [];

for fi = 1:length(files)
    %% fetching signal and annotation files
    % fi = 1;
    %fi
    % to check if annpoation file does not exsit
    if(isempty(dir(['./' alarm '/Detected_new/' files{fi},'.ann'])))
       continue;
    end
    
    [~,signal,Fs,siginfo] = rdmat(files{fi});
    %fprintf('Subject %s\n',files{fi});
    
    % to check if annpoation file is empty
    try 
        temp = dlmread(['./' alarm '/Detected_new/' files{fi},'.ann']);
    catch 
        continue;
    end
    pk_locs = temp(:,1)';
    pk_types = temp(:,2)';

    % Resample signal to 125Hz
    Fs=Fs(1);
    if Fs~=125
        signal=resample(signal,125,Fs);
        Fs=125;
    end

    xn = Fs*5*60; % alarm position
    x0 = xn-Fs*16+1; % 16s before the alarm
    signal = signal(x0:xn,:);

     description=squeeze(struct2cell(siginfo));
     description=description(4,:);
     [ecg1_orig, ecg2_orig, abp, ppg, Fs, names] = preprocess(signal, description, Fs);
     if(isempty(abp))
         abp = zeros(size(ecg1_orig));
     end
     if(isempty(ppg))
         ppg = zeros(size(ecg1_orig));
     end          

    %% separating each period of the waveform
    % Phase Wrapping

    signal_compiled = Phase_Wrapper( ecg1_orig, ecg2_orig, abp, ppg, pk_locs, Fs );
    
    signal_compiled = [signal_compiled,pk_types(2:end-1)'];

    phase_wrapped_struct(fi).signal_compiled = signal_compiled;
    phase_wrapped_struct(fi).peaks = pk_locs(2:end-1);
    phase_wrapped_struct(fi).file = files{fi};
    phase_wrapped_struct(fi).alarm_type = alarm;
    phase_wrapped_struct(fi).alarm_case = truth{fi};
    
end

%removing all entries that were empty due to annotation file problems
%(missing or empty)
phase_wrapped_struct(cellfun(@isempty,{phase_wrapped_struct.file})) = [];


%% combining
main_wrapped_signal = [];
main_fileName = [];
main_peak = [];

for ii = 1:length(phase_wrapped_struct)
    %ii
    for jj = 1:size(phase_wrapped_struct(ii).signal_compiled,1)
        main_wrapped_signal = [main_wrapped_signal;phase_wrapped_struct(ii).signal_compiled(jj,:)];
        main_fileName = [main_fileName;phase_wrapped_struct(ii).file];
        main_peak = [main_peak;phase_wrapped_struct(ii).peaks(jj)];
    end
end
wrapped_data.signal = main_wrapped_signal;
wrapped_data.fileName = main_fileName;
wrapped_data.peak = main_peak;

save([alarm '_wrapped_data_new_norm.mat'],'wrapped_data');
%%
        
        
    