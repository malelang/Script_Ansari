function peak_tester()
    data_dir = '.\';%'Y:\File Sharing Drops\Sardar\Physionet2015\Belle\7-14 Ashwin & Sardar\';

    fid=fopen([data_dir 'ALARMS'],'r');
    if(fid ~= -1)
        RECLIST=textscan(fid,'%s %s %d','Delimiter',',');
        fclose(fid);
    else
        error('Could not open ALARMS.txt for scoring. Exiting...')
    end

    RECORDS=RECLIST{1};
    ALARMS=RECLIST{2};
    TRUE_ALARM=RECLIST{3};
    N=length(RECORDS);
    
    Fs = 125;

    TP = [0 0];
    FN = [0 0];
    FP = 0;
    
    tic
    for i=1:N
        if(~strcmp(ALARMS{i},'Ventricular_Flutter_Fib') | ~(TRUE_ALARM(i)==1))
            continue;
        end            
        
        fname=RECORDS{i};
        disp(fname)
        
        [~,signal,Fs,siginfo]=rdmat([data_dir fname]);
        alarmResult=1;
        description=squeeze(struct2cell(siginfo));
        description=description(4,:);

        % Resample signal to 125Hz
        Fs=Fs(1);
        if Fs~=125
            signal=resample(signal,125,Fs);
            Fs=125;
        end        
        
        N_d=Fs*5*60; % alarm position
        N0_d=N_d-Fs*16+1; % 16s before the alarm

        [ecg1, ecg2, abp, ppg, Fs, names] = preprocess(signal, description, Fs);
        if(isempty(abp))
            abp = zeros(size(ecg1));
        end
        if(isempty(ppg))
            ppg = zeros(size(ecg1));
        end       
        ecg1 = ecg1(N0_d:N_d);
        ecg2 = ecg2(N0_d:N_d);
        abp = abp(N0_d:min(N_d,length(abp)));
        ppg = ppg(N0_d:min(N_d,length(ppg)));

        pks = peak_detection(ecg1, ecg2, abp, ppg, Fs);        
        
        test_file = dir([data_dir 'VF\Actual\' fname '.ann']);
        if(test_file.bytes==0)
            anns = [];
            ann_times = [];
            ann_types = [];               
        else
            anns = dlmread([data_dir 'VF\Actual\' fname '.ann']);
            ann_times = anns(:,1);
            ann_types = anns(:,2);            
        end
        
        ann_types(ann_types==0 | ann_types==1 | ann_types==2) = 1;
        ann_types(ann_types==3) = 2;
                
        temp = 0;
        detected = [];
        for j = 1:length(ann_times)
            if(min(abs(ann_times(j)-pks))>0.15*Fs)
                FN(ann_types(j)) = FN(ann_types(j)) + 1;
                if(ann_types(j)==2)
                    temp = temp+1;
                end
            end
        end
        
        for j = 1:length(pks)
            if(isempty(ann_times))
                FP = FP + 1;
            elseif(min(abs(ann_times-pks(j)))>0.15*Fs)
                FP = FP + 1;
            else
                [~,idx] = min(abs(ann_times-pks(j)));
                TP(ann_types(idx)) = TP(ann_types(idx)) + 1;
                detected = [detected ; pks(j),ann_types(idx)];
            end
        end    
        
        dlmwrite([data_dir 'VF\Detected\' fname '.ann'],detected);
%         if(temp>1)
%             asdf=1;
%         end
    end
    
    toc
    
    TP
    FP
    FN
    
    Se = 100*sum(TP)/(sum(TP)+sum(FN))
    PP = 100*sum(TP)/(sum(TP)+FP)

    Se1 = 100*TP(1)/(TP(1)+FN(1))
    PP1 = 100*TP(1)/(TP(1)+FP)

    Se2 = 100*TP(2)/(TP(2)+FN(2))
    PP2 = 100*TP(2)/(TP(2)+FP)
    
end

