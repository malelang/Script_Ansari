function signal_compiled = Phase_Wrapper( ecg1, ecg2, abp, ppg, pk_locs,Fs )
    
    %%%%%%%%%%%%%%% SE QUITA ESTA LINEA CUANDO SE HACE SIN FILTRO 
    ecg1 = ecg1 - medfilt1(ecg1,Fs);
    [~,init_pks] = findpeaks(abs(ecg1),'MinPeakDistance',Fs,'MinPeakHeight',prctile(abs(ecg1),70));
    if(median(abs(ecg1(init_pks)))>0)
        ecg1 = ecg1/median(abs(ecg1(init_pks)));
    end

    %%%%%%%%%%%%%%% SE QUITA ESTA LINEA CUANDO SE HACE SIN FILTRO 
    ecg2 = ecg2 - medfilt1(ecg2,Fs);
    [~,init_pks] = findpeaks(abs(ecg2),'MinPeakDistance',Fs,'MinPeakHeight',prctile(abs(ecg2),70));
    if(median(abs(ecg2(init_pks)))>0)
        ecg2 = ecg2/median(abs(ecg2(init_pks)));    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db8 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(abp,3,'db8'); 
%         A3 = wrcoef('a',C,L,'db8',3); % mejor linea base
%         abp = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(abp,4,'db6'); 
%         A3 = wrcoef('a',C,L,'db6',4); % mejor linea base
%         abp = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db10 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(abp,4,'db10'); 
%         A3 = wrcoef('a',C,L,'db10',4); % mejor linea base
%         abp = detrend(A3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym4 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(abp,3,'sym4'); 
%         A3 = wrcoef('a',C,L,'sym4',3); % mejor linea base
%         abp = detrend(A3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db10 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(abp,3,'sym6'); 
%         A3 = wrcoef('a',C,L,'sym6',3); % mejor linea base
%         abp = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%Savitzky-Golay Smoothing Filter%%%%%
%          abp2=sgolayfilt(abp,3,41);
%          abp = detrend(abp2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mov median Filter%%%%%        
%             A3  = movmedian(abp,11);
%             abp = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%% wavelets with thresholding %%%%%%%%%
%            [C,L] = wavedec(abp,5,'db6'); 
%            [thr,sorh,keepapp]=ddencmp('den','wv',abp);
%            A3=wdencmp('gbl',C,L,'db6',5,thr,sorh,keepapp);
%            abp = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%% EMD  %%%%%%%%%%%%%%%%%%%%%
  cleanedSignal = emd_dfadenoising (abp);
  cleanedSignal = cleanedSignal';
  abp=detrend(cleanedSignal); 

    %%%%%%%%%%%%%%% SE QUITA ESTA LINEA CUANDO SE HACE SIN FILTRO   
    %abp = abp - medfilt1(abp,Fs);
    [~,init_pks] = findpeaks(abs(abp),'MinPeakDistance',Fs,'MinPeakHeight',prctile(abp,70));
    if(median(abs(abp(init_pks)))>0)
        abp = abp/median(abs(abp(init_pks)));
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db8 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg,3,'db8'); 
%         A3 = wrcoef('a',C,L,'db8',3); % mejor linea base
%         ppg = detrend(A3);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg,4,'db6'); 
%         A3 = wrcoef('a',C,L,'db6',4); % mejor linea base
%         ppg = detrend(A3);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db10 level 4%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg,4,'db10'); 
%         A3 = wrcoef('a',C,L,'db10',4); % mejor linea base
%         ppg = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym4 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg,3,'sym4'); 
%         A3 = wrcoef('a',C,L,'sym4',3); % mejor linea base
%         ppg = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym6 level 3%%%%%%%%%%%%%
%         [C,L] = wavedec(ppg,3,'sym6'); 
%         A3 = wrcoef('a',C,L,'sym6',3); % mejor linea base
%         ppg = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%Savitzky-Golay Smoothing Filter%%%%%
%          ppg2=sgolayfilt(ppg,3,41);
%          ppg = detrend(ppg2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%Mov median Filter%%%%%        
%         A3  = movmedian(ppg,11);
%         ppg = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%% wavelets with thresholding %%%%%%%%%
%            [C,L] = wavedec(ppg,5,'db6'); 
%            [thr,sorh,keepapp]=ddencmp('den','wv',ppg);
%            A3=wdencmp('gbl',C,L,'db6',5,thr,sorh,keepapp);
%            ppg = detrend(A3);
%%%%%%%%%%%%%%%%%%%%%%%%%%% EMD  %%%%%%%%%%%%%%%%%%%%%
  cleanedSignal = emd_dfadenoising (ppg);
  cleanedSignal = cleanedSignal';
  ppg=detrend(cleanedSignal); 
           
    %%%%%%%%%%%%%%% SE QUITA ESTA LINEA CUANDO SE HACE SIN FILTRO      
    %ppg = ppg - medfilt1(ppg,Fs);
    [~,init_pks] = findpeaks(abs(ppg),'MinPeakDistance',Fs,'MinPeakHeight',prctile(ppg,70));
    if(median(abs(ppg(init_pks)))>0)
        ppg = ppg/median(abs(ppg(init_pks)));
    end    

    signal = [ecg1, ecg2, abp, ppg];

    %finding mid peaks for being of wave period
    pk_mids = round((pk_locs(1:end-1)+pk_locs(2:end))/2);
%     figure; plot(signal(:,1)); hold on; plot(pk_mids,signal(pk_mids,1),'*r');

    wave_interp = [];
    signal_compiled = [];
    for ii = 1:length(pk_mids)-1
        % check to see if wave period is lower than 30hz, minimize the
        % distance considered to form a wave period
        if(pk_locs(ii+1)-pk_mids(ii) > 125) % if left side of peak is larger than 1 sec
            pk_mids(ii);
            pk_locs(ii+1);
            widt = pk_locs(ii+1)-pk_mids(ii);
            %display('Wide left');
            pk_mids(ii) = pk_locs(ii+1)-125;
%             pause;
        elseif(pk_mids(ii+1)-pk_locs(ii+1) > 125)  % if right side of peak is larger than 1 sec
            pk_mids(ii+1);
            pk_locs(ii+1);
            widt = pk_mids(ii+1)-pk_locs(ii+1);
            %display('Wide right');
            pk_mids(ii+1) = pk_locs(ii+1)+125;
%             pause;
        end
        
        wave = signal(pk_mids(ii):pk_mids(ii+1),:);
        
        % zero padding ont he smaller side
        left_width = pk_locs(ii+1)-pk_mids(ii);
        right_width = pk_mids(ii+1)-pk_locs(ii+1);
        if(left_width < right_width)
            %display('left');
            padding = zeros(right_width - left_width,4);
            wave = [padding;wave];
        elseif(right_width < left_width)
            %display('right');
            padding = zeros(left_width - right_width,4);
            wave = [wave;padding];
        end
        
        wave_compiled = [];
        for jj = 1:size(signal,2)
            if(jj<3)
                if(signal(pk_locs(ii+1),jj)-median(wave(:,jj))<0)
                    wave(:,jj) = -wave(:,jj);
                end
            end
            wave_interp(:,jj) = interp1(linspace(1,125,length(wave(:,jj))),wave(:,jj),1:125); % interpolate to 125Hz (1 Sec)
            wave_compiled = [wave_compiled,wave_interp(:,jj)'];
        end
        signal_compiled(ii,:) = wave_compiled; % adding peak type (True/False)
    %     figure; plot(wave(:,1)); pause; close;
    end        
end

