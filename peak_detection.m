function pk_locs = peak_detection( ecg1_orig, ecg2_orig, abp_orig, ppg_orig, Fs )

     [bl_ecg,al_ecg] = butter(2,[0.5,40]*2/Fs);
%      [bl_abp,al_abp] = butter(2,[0.5,10]*2/Fs);
    ecg1 = filtfilt(bl_ecg,al_ecg,ecg1_orig'-ecg1_orig(1))+ecg1_orig(1);
    ecg2 = filtfilt(bl_ecg,al_ecg,ecg2_orig'-ecg2_orig(1))+ecg2_orig(1);

     %%%%%%%REMOVING FILTERS
%     ecg1=ecg1_orig';
%     ecg2=ecg2_orig';
    
    if(isempty(abp_orig))
        abp = zeros(size(ecg1));
    else
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db8 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(abp_orig,3,'db8'); 
%         A3 = wrcoef('a',C,L,'db8',3); % mejor linea base
%         abp = detrend(A3);
%         abp=abp';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(abp_orig,4,'db6'); 
%         A3 = wrcoef('a',C,L,'db6',4); % mejor linea base
%         abp = detrend(A3);
%         abp=abp';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db10 level 4%%%%%%%%%%%%%
        [C,L] = wavedec(abp_orig,4,'db10'); 
        A3 = wrcoef('a',C,L,'db10',4); % mejor linea base
        abp = detrend(A3);
        abp=abp';
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 5%%%%%%%%%%%%%
%         [C,L] = wavedec(abp_orig,5,'db6'); 
%         A3 = wrcoef('a',C,L,'db6',5); % mejor linea base
%         abp = detrend(A3);
%         abp=abp';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym4 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(abp_orig,3,'sym4'); 
%         A3 = wrcoef('a',C,L,'sym4',3); % mejor linea base
%         abp = detrend(A3);
%         abp=abp';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym6 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(abp_orig,3,'sym6'); 
%         A3 = wrcoef('a',C,L,'sym6',3); % mejor linea base
%         abp = detrend(A3);
%         abp=abp';
        %%%%%%%%%%%%%%%%%%%%%%%%%%Savitzky-Golay Smoothing Filter%%%%%
%          abp=sgolayfilt(abp_orig,3,41);
%         cleanedSignal = detrend(abp);
%         abp=cleanedSignal';
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mov median Filter%%%%%        
%             A3  = movmedian(abp_orig,11);
%             cleanedSignal = detrend(A3);
%             abp=cleanedSignal';
%       %%%%%%%%%%%%%%%%%%%%%%%%%% wavelets with thresholding %%%%%%%%%
%            [C,L] = wavedec(abp_orig,5,'db6'); 
%            [thr,sorh,keepapp]=ddencmp('den','wv',abp_orig);
%            A3=wdencmp('gbl',C,L,'db6',5,thr,sorh,keepapp);
%            abp = detrend(A3);
%            abp=abp';
        %%%%%%%%%%%%%%%%%%%%%%%%%%% EMD  %%%%%%%%%%%%%%%%%%%%%
%               cleanedSignal = emd_dfadenoising (abp_orig);
%               abp=detrend(cleanedSignal); 
        %%%%%%%%%%%%%%% FILTRO BUTTERWORTH
        %abp = filtfilt(bl_abp,al_abp,abp_orig'-abp_orig(1))+abp_orig(1);
        %%%%%%%%%%%%%%% SEÑAL ORIGINAL SIN FILTRO
        % abp=abp_orig';
    end
    if(isempty(ppg_orig))
        ppg = zeros(size(ecg1));
    else
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db8 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg_orig,3,'db8'); 
%         A3 = wrcoef('a',C,L,'db8',3); % mejor linea base
%         ppg = detrend(A3);
%         ppg=ppg';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg_orig,4,'db6'); 
%         A3 = wrcoef('a',C,L,'db6',4); % mejor linea base
%         ppg = detrend(A3);
%         ppg=ppg';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db10 level 4%%%%%%%%%%%%%
        [C,L] = wavedec(ppg_orig,4,'db10'); 
        A3 = wrcoef('a',C,L,'db10',4); % mejor linea base
        ppg = detrend(A3);
        ppg=ppg';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 5%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg_orig,5,'db6'); 
%         A3 = wrcoef('a',C,L,'db6',5); % mejor linea base
%         ppg = detrend(A3);
%         ppg=ppg';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym4 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg_orig,3,'sym4'); 
%         A3 = wrcoef('a',C,L,'sym4',3); % mejor linea base
%         ppg = detrend(A3);
%         ppg=ppg';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym6 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg_orig,3,'sym6'); 
%         A3 = wrcoef('a',C,L,'sym6',3); % mejor linea base
%         ppg = detrend(A3);
%         ppg=ppg';
        %%%%%%%%%%%%%%%%%%%%%%%%%%Savitzky-Golay Smoothing Filter%%%%%
%          ppg=sgolayfilt(ppg_orig,3,41);
%          cleanedSignal = detrend(ppg);
%          ppg=cleanedSignal';
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mov median Filter%%%%%        
%             A3  = movmedian(ppg_orig,11);
%             cleanedSignal = detrend(A3);
%             ppg=cleanedSignal';
%       %%%%%%%%%%%%%%%%%%%%%%%%%% wavelets with thresholding %%%%%%%%%
%            [C,L] = wavedec(ppg_orig,5,'db6'); 
%            [thr,sorh,keepapp]=ddencmp('den','wv',ppg_orig);
%            A3=wdencmp('gbl',C,L,'db6',5,thr,sorh,keepapp);
%            ppg = detrend(A3);
%            ppg=ppg';
        %%%%%%%%%%%%%%%%%%%%%%%%%%% EMD  %%%%%%%%%%%%%%%%%%%%%
%               cleanedSignal = emd_dfadenoising (ppg_orig);
%               ppg=detrend(cleanedSignal);
        %%%%%%%%%%%%%%% FILTRO BUTTERWORTH
       %ppg = filtfilt(bl_abp,al_abp,ppg_orig'-ppg_orig(1))+ppg_orig(1);
        %%%%%%%%%%%%%%% SEÑAL ORIGINAL SIN FILTRO
        %ppg=ppg_orig';
    end    
    
   ecg1 = ecg1 - median(ecg1);
    if((prctile(ecg1,99)-prctile(ecg1,1))>0)
        ecg1 = ecg1./(prctile(ecg1,99)-prctile(ecg1,1));
    end
   ecg2 = ecg2 - median(ecg2);
    if((prctile(ecg2,99)-prctile(ecg2,1))>0)
        ecg2 = ecg2./(prctile(ecg2,99)-prctile(ecg2,1));
    end
   %abp = abp - median(abp);
    if((prctile(abp,99)-prctile(abp,1))>0)
        abp = abp./(prctile(abp,99)-prctile(abp,1));
    end
   %ppg = ppg - median(ppg);
    if((prctile(ppg,99)-prctile(ppg,1))>0)
        ppg = ppg./(prctile(ppg,99)-prctile(ppg,1));
    end    
    
%     der_span = 3;
%     ecg1_der = abs([zeros(1,round(der_span/2)), ecg1(der_span+1:end)-ecg1(1:end-der_span), zeros(1,der_span-round(der_span/2))]);
%     ecg2_der = abs([zeros(1,round(der_span/2)), ecg2(der_span+1:end)-ecg2(1:end-der_span), zeros(1,der_span-round(der_span/2))]);
%     abp_der = ([zeros(1,round(der_span/2)), abp(der_span+1:end)-abp(1:end-der_span), zeros(1,der_span-round(der_span/2))]);
%     ppg_der = ([zeros(1,round(der_span/2)), ppg(der_span+1:end)-ppg(1:end-der_span), zeros(1,der_span-round(der_span/2))]);
% 
%     [~,ecg1_pks1] = findpeaks(ecg1_der,'MinPeakDistance',round(Fs/4));
%     [~,ecg2_pks1] = findpeaks(ecg2_der,'MinPeakDistance',round(Fs/4));
%     [~,abp_pks1] = findpeaks(abp_der,'MinPeakDistance',round(Fs/4));
%     abp_pks1 = abp_peak_search_neighborhood(abp,Fs,abp_pks1,0.2);
%     [~,ppg_pks1] = findpeaks(ppg_der,'MinPeakDistance',round(Fs/4));   

    ecg1_pks1 = my_ecg_peak_detection(ecg1,Fs);    
    ecg2_pks1 = my_ecg_peak_detection(ecg2,Fs);        
    abp_pks1 = my_abp_peak_detection(abp,Fs);
    ppg_pks1 = my_abp_peak_detection(ppg,Fs);
    
    ecg1_pks2 = sqrs_peak_detection(ecg1,Fs);
    ecg1_pks2 = ecg_peak_search_neighborhood(ecg1,Fs,ecg1_pks2);
    ecg2_pks2 = sqrs_peak_detection(ecg2,Fs);
    ecg2_pks2 = ecg_peak_search_neighborhood(ecg2,Fs,ecg2_pks2);
    
    ecg1_pks3 = peakdetect(ecg1,Fs)';
    ecg1_pks3 = ecg_peak_search_neighborhood(ecg1,Fs,ecg1_pks3);
    ecg2_pks3 = peakdetect(ecg2,Fs)';
    ecg2_pks3 = ecg_peak_search_neighborhood(ecg2,Fs,ecg2_pks3);

    bar_span = round(Fs*0.1);
    
    ecg1_indic1 = zeros(size(ecg1));
    ecg1_indic1(ecg1_pks1) = 1; %abs(ecg1_der(ecg1_pks1));
    ecg1_indic1 = conv(ecg1_indic1,ones(bar_span,1),'same');
    
    ecg1_indic2 = zeros(size(ecg1));
    ecg1_indic2(ecg1_pks2) = 1; %abs(ecg1_der(ecg1_pks1));
    ecg1_indic2 = conv(ecg1_indic2,ones(bar_span,1),'same');    
    
    ecg1_indic3 = zeros(size(ecg1));
    ecg1_indic3(ecg1_pks3) = 1; %abs(ecg1_der(ecg1_pks1));
    ecg1_indic3 = conv(ecg1_indic3,ones(bar_span,1),'same');    
    
    ecg2_indic1 = zeros(size(ecg2));
    ecg2_indic1(ecg2_pks1) = 1; %abs(ecg2_der(ecg2_pks1));
    ecg2_indic1 = conv(ecg2_indic1,ones(bar_span,1),'same');
    
    ecg2_indic2 = zeros(size(ecg2));
    ecg2_indic2(ecg2_pks2) = 1; %abs(ecg2_der(ecg2_pks1));
    ecg2_indic2 = conv(ecg2_indic2,ones(bar_span,1),'same');    
    
    ecg2_indic3 = zeros(size(ecg2));
    ecg2_indic3(ecg2_pks3) = 1; %abs(ecg2_der(ecg2_pks1));
    ecg2_indic3 = conv(ecg2_indic3,ones(bar_span,1),'same');
    
    ecg_indic = ecg1_indic1 + ecg1_indic2 + ecg1_indic3 + ecg2_indic1 + ecg2_indic2 + ecg2_indic3;
    
    if(isempty(abp))
        abp_pks2 = [];
    else
        abp_pks2 = abp_peak_detection(abp_orig);
        abp_pks2_q = abp_peak_quality(abp_orig,abp_pks2);
        abp_pks2_q = 1 - sum(abp_pks2_q(:,2:end),2)/9;
        abp_pks2 = abp_peak_search_neighborhood(abp',Fs,abp_pks2,0.2);
    end
    
    if(isempty(ppg))
        ppg_pks2 = [];
    else
        ppg_pks2 = ppg_peak_detection(ppg');
        ppg_pks2 = abp_peak_search_neighborhood(ppg',Fs,ppg_pks2,0.2);
    end    
    
    abp_indic1 = zeros(size(abp));
    abp_indic1(abp_pks1) = 3/2;%abs(abp_der(abp_pks));
    abp_indic1 = conv(abp_indic1,ones(bar_span,1),'same');
    
    abp_indic2 = zeros(size(abp));
    abp_indic2(abp_pks2) = abp_pks2_q*3/2;%abs(abp_der(abp_pks));
    abp_indic2 = conv(abp_indic2,ones(bar_span,1),'same'); 
    
    abp_indic = abp_indic1+abp_indic2;
    
    ppg_indic1 = zeros(size(ppg));
    ppg_indic1(ppg_pks1) = 3/2;%abs(ppg_der(ppg_pks));
    ppg_indic1 = conv(ppg_indic1,ones(bar_span,1),'same');
    
    ppg_indic2 = zeros(size(ppg));
    ppg_indic2(ppg_pks2) = 3/2;%abs(ppg_der(ppg_pks));
    ppg_indic2 = conv(ppg_indic2,ones(bar_span,1),'same');  
    
    ppg_indic = ppg_indic1 + ppg_indic2;
    
    %ecg_indic = ecg1_indic1 + ecg2_indic1;
    
    [r,lags] = xcorr(abp_indic,ecg_indic);
    [~,abp_lag] = max(r(lags>=0 & lags<Fs/2));
    abp_lag = abp_lag - 1;
    abp_indic_shifted = [abp_indic(abp_lag+1:end), zeros(1,abp_lag)];
    abp = [abp(abp_lag+1:end), zeros(1,abp_lag)];
    
    [r,lags] = xcorr(ppg_indic,ecg_indic);
    [~,ppg_lag] = max(r(lags>=0 & lags<2*Fs));
    ppg_lag = ppg_lag - 1;
    ppg_indic_shifted = [ppg_indic(ppg_lag+1:end), zeros(1,ppg_lag)];
    ppg = [ppg(ppg_lag+1:end), zeros(1,ppg_lag)];
    
    indic = ecg_indic + abp_indic_shifted + ppg_indic_shifted;

    indic = smooth(indic, bar_span);
    
    [~,pk_locs] = findpeaks(indic,'MinPeakDistance',round(Fs/4),'MinPeakHeight',1.75,'MinPeakProminence',1.75);
    
    adsf=1;
    
%     h = figure('units','normalized','outerposition',[0 0 1 1]);
%     h1 = subplot(3,1,1);
%     plot(ecg1);
%     hold all; 
% %    plot(pk_locs,ecg1(pk_locs),'o');    
%     plot(ecg2);
% %    plot(pk_locs,ecg2(pk_locs),'o');
%     pause(0.1);
%     ylm = get(h1,'ylim');
%     line([pk_locs';pk_locs'],repmat(ylm',1,length(pk_locs)),'linestyle','-','color','k');    
%     ylim(ylm);
%     h2 = subplot(3,1,2);
%     plot(abp);
%     hold all; 
%     %plot(pk_locs,abp(pk_locs),'o');    
%     plot(ppg);   
%     %plot(pk_locs,ppg(pk_locs),'o');
%     pause(0.1);
%     ylm = get(h2,'ylim');
%     line([pk_locs';pk_locs'],repmat(ylm',1,length(pk_locs)),'linestyle','-','color','k');
%     ylim(ylm);
%     h3 = subplot(3,1,3);
%     plot(indic);
%     hold all;plot(pk_locs,indic(pk_locs),'o');
%     linkaxes([h1,h2,h3],'x');        
%     
%     close(h);
end

function pk_locs = my_ecg_peak_detection(sig, Fs)
    der_span = 3;

    der = abs([zeros(1,round(der_span/2)), sig(der_span+1:end)-sig(1:end-der_span), zeros(1,der_span-round(der_span/2))]);
            
    [~,init_pks, init_trs] = myfindpeaks(der','MinPeakDistance',Fs*2);

    if(isempty(init_pks))
        pk_locs = [];
        return
    end    
    
    pks = der(init_pks);
    trs = [der(init_trs(1)), der(init_trs), der(init_trs(end))];
    
    depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]))/3;
    
    [~,pk_locs] = findpeaks(der','MinPeakDistance',round(Fs/4),'MinPeakProminence',depth_th);
    
%     der0 = [der(1:end-1).*der(2:end)<0, 0];
%     for i = 1:length(pk_locs)
%         [~,idx] = find(der0(pk_locs(i):end)==1,1,'first');
%         if(~isempty(idx))
%             pk_locs(i) = pk_locs(i)+idx-1;
%         end
%     end

end

function pk_locs = my_abp_peak_detection(sig, Fs)
    der_span = 3;

    der = [zeros(1,round(der_span/2)), sig(der_span+1:end)-sig(1:end-der_span), zeros(1,der_span-round(der_span/2))];
            
    [~,init_pks, init_trs] = myfindpeaks(der','MinPeakDistance',Fs*2);
    
    if(isempty(init_pks))
        pk_locs = [];
        return
    end
    
    pks = der(init_pks);
    trs = [der(init_trs(1)), der(init_trs), der(init_trs(end))];
    
    depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]))/3;
    
    [~,pk_locs] = findpeaks(der','MinPeakDistance',round(Fs/4),'MinPeakProminence',depth_th);
    
    der0 = [der(1:end-1)>=0 & der(2:end)<0, 0];
    for i = 1:length(pk_locs)
        [~,idx] = find(der0(pk_locs(i):end)==1,1,'first');
        if(~isempty(idx))
            pk_locs(i) = pk_locs(i)+idx-1;
        end
    end

end

function pk_locs = sqrs_peak_detection(sig,Fs)
    t = getCurrentTask(); 
    if(isempty(t))
        id = '1';
    else
        id = num2str(t.ID);  
    end

    sig = [zeros(1,1000), sig];

    count = 0;
    while(count < 5)
        try
            out = evalc(['mat2wfdb(sig'',''ecg' id ''',Fs,[],''v'')']);
            count = 5;
        catch
            count = count + 1;
            warning('mat2wfdb error.');
        end
    end
    sqrs(['ecg' id]);
    [pk_locs,type,subtype,chan,num,comments]=rdann(['ecg' id],'qrs');
            
    pk_locs = pk_locs - 1000;
    pk_locs(pk_locs<1) = [];    
    
%     if(any(type~='N' & type~='I'))
%         error('peak type');
%     end
end

function [pk_locs, beat_q] = abp_peak_detection(sig)
    sig = [zeros(1000,1) ; sig];
    
    pk_locs = wabp(sig,0,1);
    
    pk_locs = pk_locs - 1000;
    pk_locs(pk_locs<1) = [];
end

function beat_q = abp_peak_quality(sig,pk_locs)
    if length(pk_locs)>=3
        features = abpfeature(sig,pk_locs);
        beat_q = jSQI(features, pk_locs, sig);
    else
        beat_q = ones(size(pk_locs));
    end
end

function pk_locs = ppg_peak_detection(sig)
    sig = [zeros(1000,1) ; sig];
    
    pk_locs = wabp(sig,100*(1/20),1/20);
    
    pk_locs = pk_locs - 1000;
    pk_locs(pk_locs<1) = [];
end

function pk_locs = ecg_peak_search_neighborhood(sig,Fs,pks)

    if(length(pks)<2)
        pk_locs = pks;
        return;
    end

    der_span = 3;
    
    win = 1*Fs;

    der = [zeros(1,round(der_span/2)), sig(der_span+1:end)-sig(1:end-der_span), zeros(1,der_span-round(der_span/2))];
 
    pk_mids = round([1 mean([pks(1:end-1)';pks(2:end)']) length(sig)]);
    
    pk_locs = [];
    for i = 1:length(pks)
        [~,idx] = max(der(max(pks(i)-win,pk_mids(i)):min(pks(i)+win,pk_mids(i+1))));
        pk_locs = [pk_locs, idx+max(pks(i)-win,pk_mids(i))-1];
    end

end

function pk_locs = abp_peak_search_neighborhood(sig,Fs,pks,win_sec)

    pk_locs = [];
    win_len = round(win_sec*Fs);
    for i = 1:length(pks)
        [~,idx] = max(sig(pks(i):min(pks(i)+win_len,length(sig))));
        pk_locs = [pk_locs, idx+pks(i)-1];
    end

end
