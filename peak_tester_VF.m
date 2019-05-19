function feats = peak_tester_VF()
    data_dir = '.\';%'Y:\File Sharing Drops\Sardar\Physionet2015\Belle\7-14 Ashwin & Sardar\';

    fid=fopen([data_dir 'ALARMS'],'r');
    if(fid ~= -1)
        RECLIST=textscan(fid,'%s %s %d','Delimiter',',');
        fclose(fid);
    else
        error('Could not open ALARMS.txt for scoring. Exiting...')
    end
    
    alarm = 'VF';

    if(strcmp(alarm,'VF'))
        name = 'Ventricular_Flutter_Fib';
        max_freq = 8;
    elseif(strcmp(alarm,'VT'))
        name = 'Ventricular_Tachycardia';
        max_freq = 6;
    end    

    RECORDS=RECLIST{1};
    ALARMS=RECLIST{2};
    TRUE_ALARM=RECLIST{3};
    N=length(RECORDS);
    
    Fs = 125;

    TP1 = zeros(N,1);
    TP2 = zeros(N,1);
    FN1 = zeros(N,1);
    FN2 = zeros(N,1);
    FP = zeros(N,1);
    
    tic
    rnd = 1:N;%randperm(N);
    feats = [];
    parfor i=rnd
        if(~strcmp(ALARMS{i},name))% | ~TRUE_ALARM(i)==1)
            continue;
        end            
        
        fname=RECORDS{i};
%        disp(fname)
        
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

        if(strcmp(alarm,'VF'))
            %pks = peak_detection(ecg1, ecg2, abp, ppg, Fs); 
            pks = peak_detection_VF(ecg1, ecg2, abp, ppg, Fs, max_freq); 
        elseif(strcmp(alarm,'VT'))
            pks = peak_detection_VT(ecg1, ecg2, abp, ppg, Fs, max_freq); 
        end
        
%         pks = [];
% 
%         feats = [feats ; psd_feats(ecg1,Fs),psd_feats(ecg2,Fs),double(TRUE_ALARM(i))];
% 
%         if(feats(end,2)<1.5)
%             asdf=1;
%         end        
        
%         h = figure;
%         subplot(2,1,1);        
%         plot(pr1);
%         xlim([1,200]);
%         subplot(2,1,2);
%         plot(pr2);
%         xlim([1,200]);
%         
%         close(h);
        
        
%         [bl,al] = butter(2,0.5*2/Fs,'high');
%         feats = [feats ; mean(filter(bl,al,xcr1).^2),mean(filter(bl,al,xcr2).^2),double(TRUE_ALARM(i))];

        %feats = [feats ; stat_feats(ecg1',Fs), stat_feats(ecg2',Fs), double(TRUE_ALARM(i))];
%         if(feats(end,2)>0.1)
%             asdf=1;
%         end
%         h = figure;
%         subplot(2,1,1);
%         plot(periodogram(sum(ecg1_st1)-mean(sum(ecg1_st1))));
%         xlim([1 200]);
%         %plot(xcr1);
%         subplot(2,1,2);
%         plot(periodogram(sum(ecg2_st1)-mean(sum(ecg2_st1))));
%         xlim([1 200]);
%         %plot(xcr2);        
%         
%         close(h);

%         
%         h = figure;
%         subplot(3,2,1);
%         plot(ecg1);
%         subplot(3,2,2);
%         plot(ecg2);
%         subplot(3,2,3);
%         plot(sum(ecg1_st1))
%         hold all;
%         plot(sum(ecg1_st2))
%         subplot(3,2,4);
%         plot(sum(ecg2_st1))
%         hold all;
%         plot(sum(ecg2_st2))
%         subplot(3,2,5);
%         plot(abs(sum(ecg1_st2)-sum(ecg1_st1)))
%         subplot(3,2,6);
%         plot(abs(sum(ecg2_st2)-sum(ecg2_st1)))
%         
%         close(h);
        
        test_file = dir([data_dir alarm '\Actual\' fname '.ann']);
        if(test_file.bytes==0)
            anns = [];
            ann_times = [];
            ann_types = [];               
        else
            anns = dlmread([data_dir alarm '\Actual\' fname '.ann']);
            ann_times = anns(:,1);
            ann_types = anns(:,2);            
        end
        
        if(strcmp(alarm,'VF'))
            ann_types(ann_types==0 | ann_types==1 | ann_types==2) = 1;
            ann_types(ann_types==3) = 2;
        end
        
        span = Fs;
        if(isempty(ann_times))
            ann_start = 0;
            ann_end = 0;
        else
            ann_start = [1 ; round(mean([ann_times(2:end),ann_times(1:end-1)],2))];
            ann_start = max(ann_start,ann_times-span);
            ann_end = [round(mean([ann_times(2:end),ann_times(1:end-1)],2)) ; length(ecg1)];
            ann_end = min(ann_end,ann_times+span);            
        end
        
        count_FP = 0;
        count_FN = 0;
        detected = [];
        for j = 1:length(ann_times)
            if(~any(pks>=ann_start(j) & pks<ann_end(j))) %min(abs(ann_times(j)-pks))>0.15*Fs)
                if(ann_types(j)==1)
                    FN1(i) = FN1(i) + 1;
                else
                    FN2(i) = FN2(i) + 1;
                end
                if(ann_types(j)==2)
                    count_FN = count_FN+1;
                end
            else
                FP(i) = FP(i) + max(sum(pks>=ann_start(j) & pks<ann_end(j)) - 1,0);
                count_FP = count_FP+max(sum(pks>=ann_start(j) & pks<ann_end(j)) - 1,0);
                if(ann_types(j)==1)
                    TP1(i) = TP1(i) + 1;
                else
                    TP2(i) = TP2(i) + 1;
                end                 
                [~,idx] = min(abs(ann_times(j)-pks));
                detected = [detected ; pks(idx),ann_types(j)];
            end
        end
        
        for j = 1:length(pks)
            if(isempty(ann_times))
                FP(i) = FP(i) + 1;
                count_FP = count_FP+1;
            elseif(~any(pks(j)>=ann_start & pks(j)<ann_end)) 
                FP(i) = FP(i) + 1;
                count_FP = count_FP+1;
            end
%             elseif(min(abs(ann_times-pks(j)))>0.15*Fs)
%                 FP(i) = FP(i) + 1;
%             else
%                 [~,idx] = min(abs(ann_times-pks(j)));
%                 if(ann_types(idx)==1)
%                     TP1(i) = TP1(i) + 1;
%                 else
%                     TP2(i) = TP2(i) + 1;
%                 end                
%                 detected = [detected ; pks(j),ann_types(idx)];
%             end
        end    
        
        dlmwrite([data_dir alarm '\Detected_new\' fname '.ann'],detected);
%        count_FN
        if(count_FN>1)
            asdf=1;
        end
    end
    
    toc
    
    TP = [sum(TP1), sum(TP2)]
    FP = sum(FP)
    FN = [sum(FN1), sum(FN2)]
    
    Se = 100*sum(TP)/(sum(TP)+sum(FN))
    PP = 100*sum(TP)/(sum(TP)+FP)

    Se1 = 100*TP(1)/(TP(1)+FN(1))
    PP1 = 100*TP(1)/(TP(1)+FP)

    Se2 = 100*TP(2)/(TP(2)+FN(2))
    PP2 = 100*TP(2)/(TP(2)+FP)
    
    adsf=1;
    
end

function feats = feat_ext(sig,Fs)

    sig1 = medfilt1(sig,7);
    sig2 = abs(sig - sig1);
    sig1 = abs(sig1);

    feats = [];

    [~,init_pks, init_trs] = myfindpeaks(sig1,'MinPeakDistance',Fs*2);

    if(isempty(init_pks))
        feats = [feats, 0];
    else
        pks = sig1(init_pks);
        trs = [sig1(init_trs(1)), sig1(init_trs), sig1(init_trs(end))];

        depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]));

        feats = [feats, depth_th];        
    end            
    
    [~,init_pks, init_trs] = myfindpeaks(sig2,'MinPeakDistance',Fs*2);

    if(isempty(init_pks))
        feats = [feats, 0];
    else
        pks = sig2(init_pks);
        trs = [sig2(init_trs(1)), sig2(init_trs), sig2(init_trs(end))];

        depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]));

        feats = [feats, depth_th];        
    end     
    
end


function feats = stat_feats(sig,Fs)

    sig = sig - medfilt1(sig,Fs);
    
    if(std(sig)<0.00001)
        feats = nan(1,6);
        return
    end
    sig = (sig-mean(sig))/std(sig);
    
    wlen = 4*Fs;
    steplen = 1*Fs;
    kss = [];
    kurs = [];
    skews = [];
    for i = 1:steplen:length(sig)-wlen+1
        win = sig(i:i+wlen-1);
        
        win(win<prctile(win,25)-1.5*iqr(win) | win>prctile(win,75)+1.5*iqr(win)) = [];
        
        [~,p] = kstest(win);
        kss = [kss, p];
        kurs = [kurs kurtosis(win)];
        skews = [skews skewness(win)];
    end
    
    feats = [min(kss), max(kss), min(kurs), max(kurs), min(skews), max(skews)];

end
