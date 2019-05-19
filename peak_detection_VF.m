function pk_locs = peak_detection_VF( ecg1_orig, ecg2_orig, abp_orig, ppg_orig, Fs, max_freq )   

    if(std(abp_orig)<0.001)
        abp_orig = zeros(size(abp_orig));
    end

    if(std(ppg_orig)<0.001)
        ppg_orig = zeros(size(ppg_orig));
    end

     [bl_ecg,al_ecg] = butter(2,[0.5,40]*2/Fs);
%     [bl_abp,al_abp] = butter(2,[0.5,10]*2/Fs);
%     [bl_ppg,al_ppg] = butter(2,[0.5,max_freq]*2/Fs);
     ecg1 = filtfilt(bl_ecg,al_ecg,ecg1_orig'-ecg1_orig(1))+ecg1_orig(1);
     ecg2 = filtfilt(bl_ecg,al_ecg,ecg2_orig'-ecg2_orig(1))+ecg2_orig(1);
    
     %%%%%%%QUITANDOLE LOS FILTROS
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
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(abp_orig,4,'db6'); 
%         A3 = wrcoef('a',C,L,'db6',4); % mejor linea base
%         abp = detrend(A3);
%         abp=abp';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db10 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(abp_orig,4,'db10'); 
%         A3 = wrcoef('a',C,L,'db10',4); % mejor linea base
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
%       %%%%%%%%%%%%%%%%%%%%%%%%%%Savitzky-Golay Smoothing Filter%%%%%
%          abp=sgolayfilt(abp_orig,3,41);
%          cleanedSignal = detrend(abp);
%          abp=cleanedSignal';
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
              cleanedSignal = emd_dfadenoising (abp_orig);
              abp=detrend(cleanedSignal);

        %%%%%%%%%%%%%%% FILTRO BUTTERWORTH
        %abp = filtfilt(bl_abp,al_abp,abp_orig'-abp_orig(1))+abp_orig(1);
        %%%%%%%%%%%%%%% SEÑAL ORIGINAL SIN FILTRO
        %abp=abp_orig';
    end
    if(isempty(ppg_orig))
        ppg = zeros(size(ecg1));
    else
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db8 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg_orig,3,'db8'); 
%         A3 = wrcoef('a',C,L,'db8',3); % mejor linea base
%         ppg = detrend(A3);
%         ppg=ppg';
 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg_orig,4,'db8'); 
%         A3 = wrcoef('a',C,L,'db6',4); % mejor linea base
%         ppg = detrend(A3);
%         ppg=ppg';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db10 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg_orig,4,'db10'); 
%         A3 = wrcoef('a',C,L,'db10',4); % mejor linea base
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
%       %%%%%%%%%%%%%%%%%%%%%%%%%Savitzky-Golay Smoothing Filter%%%%%
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
              cleanedSignal = emd_dfadenoising (ppg_orig);
              ppg=detrend(cleanedSignal);

        %%%%%%%%%%%%%%% FILTRO BUTTERWORTH
        %ppg = filtfilt(bl_abp,al_abp,ppg_orig'-ppg_orig(1))+ppg_orig(1);
        %%%%%%%%%%%%%%% SEÑAL ORIGINAL SIN FILTRO
        %ppg=ppg_orig';
    end    
    
    %%%%%%%%%%%%%%% SE QUITA ESTA LINEA CUANDO SE HACE SIN FILTRO 
    ecg1 = ecg1 - median(ecg1);
    if((prctile(ecg1,99)-prctile(ecg1,1))>0)
        ecg1 = ecg1./(prctile(ecg1,99)-prctile(ecg1,1));
    end
    
    %%%%%%%%%%%%%%% SE QUITA ESTA LINEA CUANDO SE HACE SIN FILTRO 
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

    % Compute the ST of the ECG signals
    
    ecg1_ = ecg1.*[linspace(0,1,round(Fs/2)), ones(1,length(ecg1)-2*round(Fs/2)), linspace(1,0,round(Fs/2))];
    ecg1_st = abs(stran(ecg1_-mean(ecg1_)));    
    ecg1_st_sum = mean(ecg1_st); 
    
    ecg2_ = ecg2.*[linspace(0,1,round(Fs/2)), ones(1,length(ecg2)-2*round(Fs/2)), linspace(1,0,round(Fs/2))];
    ecg2_st = abs(stran(ecg2_-mean(ecg2_)));    
    ecg2_st_sum = mean(ecg2_st);     
    
    
    % ECG Peak Detection
    
    [ecg1_pks1,ecg1_pks1_q] = my_ecg_peak_detection(ecg1,Fs,max_freq,ecg1_st_sum);    
    [ecg2_pks1,ecg2_pks1_q] = my_ecg_peak_detection(ecg2,Fs,max_freq,ecg2_st_sum);                  
    
%     [ecg1_pks2,ecg1_pks2_q] = sqrs_peak_detection(ecg1,Fs,max_freq,ecg1_st_sum);
%     %ecg1_pks2 = ecg_peak_search_neighborhood(ecg1,Fs,ecg1_pks2);
%     [ecg2_pks2,ecg2_pks2_q] = sqrs_peak_detection(ecg2,Fs,max_freq,ecg2_st_sum);
%     %ecg2_pks2 = ecg_peak_search_neighborhood(ecg2,Fs,ecg2_pks2);
%     
%     [ecg1_pks3,ecg1_pks3_q] = sadeghi_peak_detection(ecg1,Fs,max_freq,ecg1_st_sum);
%     %ecg1_pks3 = ecg_peak_search_neighborhood(ecg1,Fs,ecg1_pks3);
%     [ecg2_pks3,ecg2_pks3_q] = sadeghi_peak_detection(ecg2,Fs,max_freq,ecg2_st_sum);
%     %ecg2_pks3 = ecg_peak_search_neighborhood(ecg2,Fs,ecg2_pks3);
    
    [ecg1_pks4,ecg1_pks4_q] = ecg_st_peak_detection(ecg1,Fs,max_freq,ecg1_st_sum);  
    %ecg1_pks4 = ecg_peak_search_neighborhood(ecg1,Fs,ecg1_pks4);
    [ecg2_pks4,ecg2_pks4_q] = ecg_st_peak_detection(ecg2,Fs,max_freq,ecg2_st_sum);    
    %ecg2_pks4 = ecg_peak_search_neighborhood(ecg2,Fs,ecg2_pks4);

    [ecg1_pks5,ecg1_pks5_q] = my_ecg_peak_detection2(ecg1,Fs,max_freq,ecg1_st_sum);  
    %ecg1_pks4 = ecg_peak_search_neighborhood(ecg1,Fs,ecg1_pks4);
    if(corr(ecg1',ecg2')<0)
        [ecg2_pks5,ecg2_pks5_q] = my_ecg_peak_detection2(-ecg2,Fs,max_freq,ecg2_st_sum);    
    else
        [ecg2_pks5,ecg2_pks5_q] = my_ecg_peak_detection2(ecg2,Fs,max_freq,ecg2_st_sum);    
    end
    %ecg2_pks4 = ecg_peak_search_neighborhood(ecg2,Fs,ecg2_pks4);
    
    
%     h = figure;
%     subplot(2,2,1);
%     plot(ecg1);
%     hold all;plot(ecg1_pks4,ecg1(ecg1_pks4),'o');
%     subplot(2,2,2);
%     plot(ecg2);
%     hold all;plot(ecg2_pks4,ecg2(ecg2_pks4),'o');
%     subplot(2,2,3);
%     plot(ecg1_pks4,ecg1_pks4_q);%abp_pks2_q);
%     ylim([0,1]);
%     subplot(2,2,4);
%     plot(ecg2_pks4,ecg2_pks4_q);
%     ylim([0,1]);
%     
%     close(h);   
        
    ecg1_indic = indicator(ecg1,ecg1_pks1,ecg1_pks1_q,Fs) + indicator(ecg1,ecg1_pks4,ecg1_pks4_q,Fs) + indicator(ecg1,ecg1_pks5,ecg1_pks5_q,Fs);
    ecg2_indic = indicator(ecg2,ecg2_pks1,ecg2_pks1_q,Fs) + indicator(ecg2,ecg2_pks4,ecg2_pks4_q,Fs) + indicator(ecg2,ecg2_pks5,ecg2_pks5_q,Fs);        

%     ecg1_indic = indicator(ecg1,ecg1_pks1,ecg1_pks1_q,Fs) + indicator(ecg1,ecg1_pks2,ecg1_pks2_q,Fs) + ...
%                  indicator(ecg1,ecg1_pks3,ecg1_pks3_q,Fs) + indicator(ecg1,ecg1_pks4,ecg1_pks4_q,Fs);
%     ecg2_indic = indicator(ecg2,ecg2_pks1,ecg2_pks1_q,Fs) + indicator(ecg2,ecg2_pks2,ecg2_pks2_q,Fs) + ...
%                  indicator(ecg2,ecg2_pks3,ecg2_pks3_q,Fs) + indicator(ecg2,ecg2_pks4,ecg2_pks4_q,Fs);
    
    ecg_indic = ecg1_indic + ecg2_indic;           
           
           
    % ABP Peak Detection
    
%     if(isempty(abp))
%         abp_pks1 = [];
%         abp_pks2 = [];
%     else
%         [abp_pks1, abp_pks1_q] = my_abp_peak_detection(abp,Fs,max_freq);
%         [abp_pks2, abp_pks2_q] = abp_peak_detection(abp_orig,Fs,max_freq);               
%     end
%     
%     abp_indic = indicator(abp,abp_pks1,abp_pks1_q,Fs) + indicator(abp,abp_pks2,abp_pks2_q,Fs);
%     
%     % PPG Peak Detection
%     
%     if(isempty(ppg))
%         ppg_pks1 = [];
%         ppg_pks2 = [];
%     else
%         [ppg_pks1, ppg_pks1_q] = my_ppg_peak_detection(ppg,Fs,max_freq);
%         [ppg_pks2, ppg_pks2_q] = ppg_peak_detection(ppg',Fs,max_freq);
%     end   
%     
%     ppg_indic = indicator(ppg,ppg_pks1,ppg_pks1_q,Fs) + indicator(ppg,ppg_pks2,ppg_pks2_q,Fs);
    
%     abp_indic1 = zeros(size(abp));
%     abp_indic1(abp_pks1) = 3/2;%abs(abp_der(abp_pks));
%     abp_indic1 = conv(abp_indic1,ones(bar_span,1),'same');
%     
%     abp_indic2 = zeros(size(abp));
%     abp_indic2(abp_pks2) = abp_pks2_q*3/2;%abs(abp_der(abp_pks));
%     abp_indic2 = conv(abp_indic2,ones(bar_span,1),'same'); 
    
%     abp_indic3 = zeros(size(abp));
%     abp_indic3(abp_pks3) = 3/2;
%     abp_indic3 = conv(abp_indic3,ones(bar_span,1),'same');     
    
%    abp_indic = abp_indic1 + abp_indic2;% + abp_indic3;
    
%     ppg_indic1 = zeros(size(ppg));
%     ppg_indic1(ppg_pks1) = 3/2;%abs(ppg_der(ppg_pks));
%     ppg_indic1 = conv(ppg_indic1,ones(bar_span,1),'same');
%     
%     ppg_indic2 = zeros(size(ppg));
%     ppg_indic2(ppg_pks2) = 3/2;%abs(ppg_der(ppg_pks));
%     ppg_indic2 = conv(ppg_indic2,ones(bar_span,1),'same');  
    
%     ppg_indic3 = zeros(size(ppg));
%     ppg_indic3(ppg_pks3) = 3/2;
%     ppg_indic3 = conv(ppg_indic3,ones(bar_span,1),'same');   
    
%     h = figure;
%     subplot(2,2,1);
%     plot(abp);
%     hold all;plot(abp_pks1,abp(abp_pks1),'o');
%     subplot(2,2,2);
%     plot(ppg);
%     hold all;plot(ppg_pks1,ppg(ppg_pks1),'o');
%     subplot(2,2,3);
%     plot(abp_pks1,abp_peak_quality(abp_orig,abp_pks1));%abp_pks2_q);
%     ylim([0,1]);
%     subplot(2,2,4);
%     plot(ppg_pks1,ppg_peak_quality(ppg,ppg_pks1,max_freq));
%     ylim([0,1]);
%     
%     close(h);    
    
%    ppg_indic = ppg_indic1 + ppg_indic2;% + ppg_indic3;
    
    %ecg_indic = ecg1_indic1 + ecg2_indic1;
    
%     [r,lags] = xcorr(abp_indic,ecg_indic);
%     [~,abp_lag] = max(r(lags>=0 & lags<Fs/2));
%     abp_lag = abp_lag - 1;
%     abp_indic_shifted = [abp_indic(abp_lag+1:end), zeros(1,abp_lag)];
%     abp = [abp(abp_lag+1:end), zeros(1,abp_lag)];
%     
%     [r,lags] = xcorr(ppg_indic,ecg_indic);
%     [~,ppg_lag] = max(r(lags>=0 & lags<2*Fs));
%     ppg_lag = ppg_lag - 1;
%     ppg_indic_shifted = [ppg_indic(ppg_lag+1:end), zeros(1,ppg_lag)];
%     ppg = [ppg(ppg_lag+1:end), zeros(1,ppg_lag)];
    
    indic = ecg_indic;% + abp_indic_shifted + ppg_indic_shifted;

    bar_span = round(Fs*0.1);
    
    indic = smooth(indic, bar_span);
    
    [~,pk_locs] = findpeaks(indic,'MinPeakDistance',round(Fs/max_freq),'MinPeakHeight',1.2,'MinPeakProminence',1);
    
    pk_locs = search_neighborhood(ecg1, pk_locs, floor(Fs/max_freq/2));
    
%     h = figure;
%     ah1 = subplot(4,2,1);
%     plot(ecg1);
%     hold all;plot(pk_locs,ecg1(pk_locs),'o');
%     ah2 = subplot(4,2,2);
%     plot(ecg2);
%     hold all;plot(pk_locs,ecg2(pk_locs),'o');
%     ah3 = subplot(4,2,3);
%     plot(ecg1_indic);
%     hold all;plot(pk_locs,ecg1_indic(pk_locs),'o');
%     ah4 = subplot(4,2,4);
%     plot(ecg2_indic);
%     hold all;plot(pk_locs,ecg2_indic(pk_locs),'o');
% 
%     ah5 = subplot(4,2,5);
%     plot(abp);
%     hold all;plot(pk_locs,abp(pk_locs),'o');
%     ah6 = subplot(4,2,6);
%     plot(ppg);
%     hold all;plot(pk_locs,ppg(pk_locs),'o');
%     ah7 = subplot(4,2,7);
%     plot(abp_indic_shifted);
%     hold all;plot(pk_locs,abp_indic_shifted(pk_locs),'o');
%     ah8 = subplot(4,2,8);
%     plot(ppg_indic_shifted);
%     hold all;plot(pk_locs,ppg_indic_shifted(pk_locs),'o');    
%     
%     linkaxes([ah1,ah2,ah3,ah4,ah5,ah6,ah7,ah8], 'x');
%     
%     close(h); 
    
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

function indic = indicator(sig,pk_locs,pk_q,Fs)
    bar_span = round(Fs*0.1);
    
    indic = zeros(size(sig));
    indic(pk_locs) = pk_q;
    indic = conv(indic,ones(bar_span,1),'same');
end

function [pk_locs, pk_q] = my_ecg_peak_detection(sig, Fs, max_freq, st_sum)
    der_span = 3;

    der = abs([zeros(1,round(der_span/2)), sig(der_span+1:end)-sig(1:end-der_span), zeros(1,der_span-round(der_span/2))]);
            
    [~,init_pks, init_trs] = myfindpeaks(der','MinPeakDistance',Fs*2);

    if(isempty(init_pks))
        pk_locs = [];
        pk_q = [];
        return
    end    
    
    pks = der(init_pks);
    trs = [der(init_trs(1)), der(init_trs), der(init_trs(end))];
    
    depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]))/3;
    
    [~,pk_locs] = findpeaks(der','MinPeakDistance',round(Fs/max_freq),'MinPeakProminence',depth_th);
    
    pk_locs = search_neighborhood(sig, pk_locs, floor(Fs/max_freq/2));
    
    pk_q = peak_quality(sig,pk_locs,floor(Fs/max_freq/2));
    
%     der0 = [der(1:end-1).*der(2:end)<0, 0];
%     for i = 1:length(pk_locs)
%         [~,idx] = find(der0(pk_locs(i):end)==1,1,'first');
%         if(~isempty(idx))
%             pk_locs(i) = pk_locs(i)+idx-1;
%         end
%     end

end


function [pk_locs, pk_q] = my_ecg_peak_detection2(sig, Fs, max_freq, st_sum)
            
    [~,init_pks, init_trs] = myfindpeaks(sig,'MinPeakDistance',Fs*2);

    if(isempty(init_pks))
        pk_locs = [];
        pk_q = [];
        return
    end    
    
    pks = sig(init_pks);
    trs = [sig(init_trs(1)), sig(init_trs), sig(init_trs(end))];
    
    depth_th = min(mean([pks-trs(1:end-1);pks-trs(2:end)]))/4;
    
    [~,pk_locs] = findpeaks(sig,'MinPeakDistance',round(Fs/max_freq),'MinPeakProminence',depth_th);
    
%    pk_locs = search_neighborhood(sig, pk_locs, floor(Fs/max_freq/2));
    
    pk_q = peak_quality(sig,pk_locs,floor(Fs/max_freq/2));
    
end


function [pk_locs, pk_q] = sqrs_peak_detection(sig,Fs,max_freq,st_sum)
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
    
    pk_locs = search_neighborhood(sig, pk_locs, floor(Fs/max_freq/2));
    
    pk_q = peak_quality(sig,pk_locs,floor(Fs/max_freq/2));
    
%     if(any(type~='N' & type~='I'))
%         error('peak type');
%     end
end

function [pk_locs, pk_q] = sadeghi_peak_detection(sig,Fs,max_freq,st_sum)
    pk_locs = peakdetect(sig,Fs)';
    
    pk_locs = search_neighborhood(sig, pk_locs, floor(Fs/max_freq/2));
    
    pk_q = peak_quality(sig,pk_locs,floor(Fs/max_freq/2));
end

function [pk_locs, pk_q] = ecg_st_peak_detection(sig, Fs, max_freq, st_sum)

    [~,init_pks, init_trs] = myfindpeaks(st_sum,'MinPeakDistance',Fs*2);
    
    if(isempty(init_pks))
        pk_locs = [];
        pk_q = [];
        return
    elseif(length(init_pks)==1)
        [~,idx1] = min(st_sum(1:init_pks));
        [~,idx2] = min(st_sum(init_pks+1:end));
        init_trs = [idx1, idx2+init_pks];
    end
    
    pks = st_sum(init_pks);
    trs = [st_sum(init_trs(1)), st_sum(init_trs), st_sum(init_trs(end))];
    
    depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]))/5;
    
    [~,pk_locs] = findpeaks(st_sum','MinPeakDistance',round(Fs/max_freq),'MinPeakProminence',depth_th);
    
    pk_locs = search_neighborhood(sig, pk_locs, floor(Fs/max_freq/2));      
        
    pk_q = peak_quality(sig,pk_locs,floor(Fs/max_freq/2));
        
end





function [pk_locs, pk_q] = my_abp_peak_detection(sig, Fs, max_freq)
    der_span = 3;

    der = [zeros(1,round(der_span/2)), sig(der_span+1:end)-sig(1:end-der_span), zeros(1,der_span-round(der_span/2))];
            
    [~,init_pks, init_trs] = myfindpeaks(der','MinPeakDistance',Fs*2);
    
    if(isempty(init_pks))
        pk_locs = [];
        pk_q = [];
        return
    end
    
    pks = der(init_pks);
    trs = [der(init_trs(1)), der(init_trs), der(init_trs(end))];
    
    depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]))/3;
    
    [~,pk_locs] = findpeaks(der','MinPeakDistance',round(Fs/max_freq),'MinPeakProminence',depth_th);

    pk_q = abp_peak_quality(sig,pk_locs);
    
    %pk_locs = search_neighborhood(sig, pk_locs, floor(Fs/max_freq/2));  
    [pk_locs,deleted] = abp_search_neighborhood(sig,Fs,pk_locs,1,floor(Fs/max_freq/2));            
    pk_q(deleted) = [];
    
%     der0 = [der(1:end-1)>=0 & der(2:end)<0, 0];
%     for i = 1:length(pk_locs)
%         [~,idx] = find(der0(pk_locs(i):end)==1,1,'first');
%         if(~isempty(idx))
%             pk_locs(i) = pk_locs(i)+idx-1;
%         end
%     end

end

function [pk_locs, pk_q] = abp_peak_detection(sig,Fs,max_freq)
    sig2 = [zeros(1000,1) ; sig];
    
    pk_locs = wabp(sig2,0,1);
    
    pk_locs = pk_locs - 1000;
    pk_locs(pk_locs<1) = [];                
    
    pk_q = abp_peak_quality(sig',pk_locs);
    
    %pk_locs = search_neighborhood(sig', pk_locs,floor(Fs/max_freq/2));    
    [pk_locs,deleted] = abp_search_neighborhood(sig',Fs,pk_locs,1,floor(Fs/max_freq/2));
    pk_q(deleted) = [];
end

function beat_q = abp_peak_quality(sig,pk_locs)
    if length(pk_locs)>=3
        features = abpfeature(sig,pk_locs);
        beat_q = jSQI(features, pk_locs, sig);
        beat_q = 1 - sum(beat_q(:,2:end),2)/9;
        beat_q = [beat_q(1:size(features,1))', repmat(beat_q(size(features,1)),1,length(pk_locs)-size(features,1))];
    else
        beat_q = ones(size(pk_locs));
    end
end

function pk_q = abp_peak_quality2(sig,pk_locs,Fs,max_freq)

    span = floor(Fs/max_freq/2); 

    if(isempty(pk_locs))
        pk_q = [];
        return
    end

    pk_q = zeros(size(pk_locs));        
    for i = 2:length(pk_locs)-1
        pr1 = sig(max(pk_locs(i-1)-span,1):min(pk_locs(i-1)+span,length(sig)));
        pr1 = [repmat(pr1(1),1,span-pk_locs(i-1)), pr1, repmat(pr1(end),1,span-(length(sig)-pk_locs(i-1)))];

        pr2 = sig(max(pk_locs(i)-span,1):min(pk_locs(i)+span,length(sig)));
        pr2 = [repmat(pr2(1),1,span-pk_locs(i)), pr2, repmat(pr2(end),1,span-(length(sig)-pk_locs(i)))];
        
        pr3 = sig(max(pk_locs(i+1)-span,1):min(pk_locs(i+1)+span,length(sig)));
        pr3 = [repmat(pr3(1),1,span-pk_locs(i+1)), pr3, repmat(pr3(end),1,span-(length(sig)-pk_locs(i+1)))];        
        
        pk_q(i) = max([abs(corr(pr1',pr2')), abs(corr(pr2',pr3'))]);
    end

    pr2 = sig(max(pk_locs(1)-span,1):min(pk_locs(1)+span,length(sig)));
    pr2 = [repmat(pr2(1),1,span-pk_locs(1)), pr2, repmat(pr2(end),1,span-(length(sig)-pk_locs(1)))];

    pr3 = sig(max(pk_locs(1+1)-span,1):min(pk_locs(1+1)+span,length(sig)));
    pr3 = [repmat(pr3(1),1,span-pk_locs(1+1)), pr3, repmat(pr3(end),1,span-(length(sig)-pk_locs(1+1)))];         
    
    pk_q(1) = abs(corr(pr2',pr3'));

    pr1 = sig(max(pk_locs(end-1)-span,1):min(pk_locs(end-1)+span,length(sig)));
    pr1 = [repmat(pr1(1),1,span-pk_locs(end-1)), pr1, repmat(pr1(end),1,span-(length(sig)-pk_locs(end-1)))];

    pr2 = sig(max(pk_locs(end)-span,1):min(pk_locs(end)+span,length(sig)));
    pr2 = [repmat(pr2(1),1,span-pk_locs(end)), pr2, repmat(pr2(end),1,span-(length(sig)-pk_locs(end)))];
    
    pk_q(end) = abs(corr(pr1',pr2'));
end

function [pk_locs, pk_q] = my_ppg_peak_detection(sig, Fs, max_freq)
    der_span = 3;

    der = [zeros(1,round(der_span/2)), sig(der_span+1:end)-sig(1:end-der_span), zeros(1,der_span-round(der_span/2))];
            
    [~,init_pks, init_trs] = myfindpeaks(der','MinPeakDistance',Fs*2);
    
    if(isempty(init_pks))
        pk_locs = [];
        pk_q = [];
        return
    end
    
    pks = der(init_pks);
    trs = [der(init_trs(1)), der(init_trs), der(init_trs(end))];
    
    depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]))/3;
    
    [~,pk_locs] = findpeaks(der','MinPeakDistance',round(Fs/max_freq),'MinPeakProminence',depth_th);
    
    zr_cr = find(der(2:end)>0 & der(1:end-1)<=0);

    del = [];
    for i = 1:length(pk_locs)
        [~,idx] = find(zr_cr<pk_locs(i),1,'last');
        if(isempty(idx))
            del = [del i];
            continue;
        end
        pk_locs(i) = zr_cr(idx);
    end
    pk_locs(del) = [];        

    pk_q = ppg_peak_quality(sig,pk_locs,max_freq);
    
    %pk_locs = search_neighborhood(sig, pk_locs, floor(Fs/max_freq/2));  
    [pk_locs,deleted] = abp_search_neighborhood(sig,Fs,pk_locs,1,floor(Fs/max_freq/2));            
    pk_q(deleted) = [];
    
%     der0 = [der(1:end-1)>=0 & der(2:end)<0, 0];
%     for i = 1:length(pk_locs)
%         [~,idx] = find(der0(pk_locs(i):end)==1,1,'first');
%         if(~isempty(idx))
%             pk_locs(i) = pk_locs(i)+idx-1;
%         end
%     end

end

function [pk_locs,pk_q] = ppg_peak_detection(sig, Fs, max_freq)
    sig2 = [zeros(1000,1) ; sig];
    
    pk_locs = wabp(sig2,100*(1/20),1/20);
    
    pk_locs = pk_locs - 1000;
    pk_locs(pk_locs<1) = [];
    
    pk_q = ppg_peak_quality(sig',pk_locs,max_freq);        
    
    %pk_locs = search_neighborhood(sig', pk_locs,floor(Fs/max_freq/2));        
    [pk_locs,deleted] = abp_search_neighborhood(sig',Fs,pk_locs,1,floor(Fs/max_freq/2));
    pk_q(deleted) = [];
end

function beat_q = ppg_peak_quality(sig,pk_locs,max_freq)
    if(~isempty(pk_locs))
        beat_q = ppgSQI2(sig,pk_locs,max_freq);
    else
        beat_q = [];
    end
end

function pk_locs = ecg_search_neighborhood(sig,Fs,pks)

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

function [pk_locs,del] = abp_search_neighborhood(sig,Fs,pk_locs,win_sec,span)
    del = [];

    if(isempty(pk_locs))
        return;
    end
    
    win_len = round(win_sec*Fs);
%     for i = 1:length(pks)-1
%         [~,idx] = max(sig(max(pks(i)-span,1):min([pks(i)+win_len,length(sig),pks(i+1)])));
%         pk_locs = [pk_locs, idx+max(pks(i)-span,1)-1];
%     end
%     i = length(pks);
%     [~,idx] = max(sig(max(pks(i)-span,1):min([pks(i)+win_len,length(sig)])));
%     pk_locs = [pk_locs, idx+max(pks(i)-span,1)-1];    
    
    
    dr = [0 sig(2:end)-sig(1:end-1)];
    zr_cr = find(dr(2:end)<0 & dr(1:end-1)>=0);

    for i = 1:length(pk_locs)
        if(i == length(pk_locs))
            [~,idx] = find(zr_cr>max(pk_locs(i)-span,1) & zr_cr<min([pk_locs(i)+win_len,length(sig)]));
        else
            [~,idx] = find(zr_cr>max(pk_locs(i)-span,1) & zr_cr<min([pk_locs(i)+win_len,length(sig),pk_locs(i+1)-1]));
        end
        if(isempty(idx))
            del = [del i];
            continue;
        end
        [~,idx2] = max(abs(sig(zr_cr(idx))));
        pk_locs(i) = zr_cr(idx(idx2));
    end
    pk_locs(del) = [];    

end



function [pk_locs, pk_q] = abp_st_peak_detection(sig, Fs, max_freq)
    
    sig2 = sig - median(sig);
    sig2 = sig2.*[linspace(0,1,round(Fs/2)), ones(1,length(sig2)-2*round(Fs/2)), linspace(1,0,round(Fs/2))];

    sig_st = abs(st(sig2-mean(sig2)));
    
    st_sum = mean(sig_st);

    [~,init_pks, init_trs] = myfindpeaks(st_sum,'MinPeakDistance',Fs*2);
    
    if(isempty(init_pks))
        pk_locs = [];
        pk_q = [];
        return
    elseif(length(init_pks)==1)
        [~,idx1] = min(st_sum(1:init_pks));
        [~,idx2] = min(st_sum(init_pks+1:end));
        init_trs = [idx1, idx2+init_pks];
    end
    
    pks = st_sum(init_pks);
    trs = [st_sum(init_trs(1)), st_sum(init_trs), st_sum(init_trs(end))];
    
    depth_th = median(mean([pks-trs(1:end-1);pks-trs(2:end)]))/5;
    
    [~,pk_locs] = findpeaks(st_sum','MinPeakDistance',round(Fs/max_freq),'MinPeakProminence',depth_th);
    
    pk_locs = search_neighborhood(sig, pk_locs, floor(Fs/max_freq/2));
        
    pk_q = zeros(size(pk_locs));        
    for i = 2:length(pk_locs)-1
        pr1 = st_sum(max(pk_locs(i-1)-span,1):min(pk_locs(i-1)+span,length(st_sum)));
        pr1 = [repmat(pr1(1),1,span-pk_locs(i-1)), pr1, repmat(pr1(end),1,span-(length(st_sum)-pk_locs(i-1)))];

        pr2 = st_sum(max(pk_locs(i)-span,1):min(pk_locs(i)+span,length(st_sum)));
        pr2 = [repmat(pr2(1),1,span-pk_locs(i)), pr2, repmat(pr2(end),1,span-(length(st_sum)-pk_locs(i)))];
        
        pr3 = st_sum(max(pk_locs(i+1)-span,1):min(pk_locs(i+1)+span,length(st_sum)));
        pr3 = [repmat(pr3(1),1,span-pk_locs(i+1)), pr3, repmat(pr3(end),1,span-(length(st_sum)-pk_locs(i+1)))];        
        
        pk_q(i) = max([abs(corr(pr1',pr2')), abs(corr(pr2',pr3'))]);
    end

    pr2 = st_sum(max(pk_locs(1)-span,1):min(pk_locs(1)+span,length(st_sum)));
    pr2 = [repmat(pr2(1),1,span-pk_locs(1)), pr2, repmat(pr2(end),1,span-(length(st_sum)-pk_locs(1)))];

    pr3 = st_sum(max(pk_locs(1+1)-span,1):min(pk_locs(1+1)+span,length(st_sum)));
    pr3 = [repmat(pr3(1),1,span-pk_locs(1+1)), pr3, repmat(pr3(end),1,span-(length(st_sum)-pk_locs(1+1)))];         
    
    pk_q(1) = abs(corr(pr2',pr3'));

    pr1 = st_sum(max(pk_locs(end-1)-span,1):min(pk_locs(end-1)+span,length(st_sum)));
    pr1 = [repmat(pr1(1),1,span-pk_locs(end-1)), pr1, repmat(pr1(end),1,span-(length(st_sum)-pk_locs(end-1)))];

    pr2 = st_sum(max(pk_locs(end)-span,1):min(pk_locs(end)+span,length(st_sum)));
    pr2 = [repmat(pr2(1),1,span-pk_locs(end)), pr2, repmat(pr2(end),1,span-(length(st_sum)-pk_locs(end)))];
    
    pk_q(end) = abs(corr(pr1',pr2'));
        
end


function pk_locs = search_neighborhood(sig, pk_locs, span)
    dr = [0 sig(2:end)-sig(1:end-1)];
    zr_cr = find((dr(2:end).*dr(1:end-1))<0);

    del = [];
    for i = 1:length(pk_locs)
        [~,idx] = find(abs(zr_cr-pk_locs(i))<span);
        if(isempty(idx))
            %del = [del i];
            continue;
        end
        [~,idx2] = max(abs(sig(zr_cr(idx))));
        pk_locs(i) = zr_cr(idx(idx2));
    end
    pk_locs(del) = [];  
end


function pk_q = peak_quality(sig,pk_locs,span)
    if(length(pk_locs)<3)
        pk_q = 0.4*ones(size(pk_locs));
        return
    end

    pk_q = zeros(size(pk_locs));        
    for i = 2:length(pk_locs)-1
        pr1 = sig(max(pk_locs(i-1)-span,1):min(pk_locs(i-1)+span,length(sig)));
        pr1 = [repmat(pr1(1),1,span-pk_locs(i-1)+1), pr1, repmat(pr1(end),1,span-(length(sig)-pk_locs(i-1)))];

        pr2 = sig(max(pk_locs(i)-span,1):min(pk_locs(i)+span,length(sig)));
        pr2 = [repmat(pr2(1),1,span-pk_locs(i)+1), pr2, repmat(pr2(end),1,span-(length(sig)-pk_locs(i)))];
        
        pr3 = sig(max(pk_locs(i+1)-span,1):min(pk_locs(i+1)+span,length(sig)));
        pr3 = [repmat(pr3(1),1,span-pk_locs(i+1)+1), pr3, repmat(pr3(end),1,span-(length(sig)-pk_locs(i+1)))];        
        
        pk_q(i) = max([abs(corr(pr1',pr2')), abs(corr(pr2',pr3'))]);
    end

    pr2 = sig(max(pk_locs(1)-span,1):min(pk_locs(1)+span,length(sig)));
    pr2 = [repmat(pr2(1),1,span-pk_locs(1)+1), pr2, repmat(pr2(end),1,span-(length(sig)-pk_locs(1)))];

    pr3 = sig(max(pk_locs(1+1)-span,1):min(pk_locs(1+1)+span,length(sig)));
    pr3 = [repmat(pr3(1),1,span-pk_locs(1+1)+1), pr3, repmat(pr3(end),1,span-(length(sig)-pk_locs(1+1)))];         
    
    pk_q(1) = abs(corr(pr2',pr3'));

    pr1 = sig(max(pk_locs(end-1)-span,1):min(pk_locs(end-1)+span,length(sig)));
    pr1 = [repmat(pr1(1),1,span-pk_locs(end-1)+1), pr1, repmat(pr1(end),1,span-(length(sig)-pk_locs(end-1)))];

    pr2 = sig(max(pk_locs(end)-span,1):min(pk_locs(end)+span,length(sig)));
    pr2 = [repmat(pr2(1),1,span-pk_locs(end)+1), pr2, repmat(pr2(end),1,span-(length(sig)-pk_locs(end)))];
    
    pk_q(end) = abs(corr(pr1',pr2'));
end