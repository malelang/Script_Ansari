function [ecg1, ecg2, abp, ppg, Fs, names] = preprocess(sig, description, Fs,numRecord)

    names = {};

    % Resample sig to 125Hz
    Fs=Fs(1);
    if Fs~=125
        sig=resample(sig,125,Fs);
        Fs=125;
    end    

    ecg_names = {'II','V','III','aVF','I','MCL','aVL','aVR'};
    ecg_indices = {};
    for i = 1:length(ecg_names)
        ind = find(strcmp(description,ecg_names{i}));
        if(length(ind)==1)
            ecg_indices = [ecg_indices, {ecg_names{i};ind}];
        end
    end    

    if(size(ecg_indices,2)==0)
        warning(['Cannot find two ECG leads. recordName = ' recordName]);
    end

    ecg1_ind = [];
    ecg2_ind = [];
    if(any(strcmp(ecg_indices(1,:),'II')))
        ecg1_ind = ecg_indices{2,strcmp(ecg_indices(1,:),'II')};
        ecg_indices(:,strcmp(ecg_indices(1,:),'II')) = [];
    end    
    if(any(strcmp(ecg_indices(1,:),'V')))
        ecg2_ind = ecg_indices{2,strcmp(ecg_indices(1,:),'V')};
        ecg_indices(:,strcmp(ecg_indices(1,:),'V')) = [];
    end
    if(isempty(ecg1_ind))
        ecg1_ind = ecg_indices{2,1};
        ecg_indices(:,1) = [];
    end
    if(isempty(ecg2_ind))
        ecg2_ind = ecg_indices{2,1};
        ecg_indices(:,1) = [];
    end
    ecg1 = sig(:,ecg1_ind);
    ecg2 = sig(:,ecg2_ind);
    
    names = {description{ecg1_ind}, description{ecg2_ind}, 'ABP', 'PLETH'};

    ecg1(isnan(ecg1)) = 0;
%     if(length(unique(ecg1))<1000 || sum(ecg1==0)>30*Fs)
%         ecg1 = [];
%     end        
    
    ecg2(isnan(ecg2)) = 0;
%     if(length(unique(ecg2))<1000 || sum(ecg2==0)>30*Fs)
%         ecg2 = [];
%     end    
    
    abp_ind = find(strcmp(description,'ABP'));    
    if(isempty(abp_ind))
        abp = [];
    else
        abp = sig(:,abp_ind);
        abp(isnan(abp)) = 0;
        signal=abp;
    end
    
    ppg_ind = find(strcmp(description,'PLETH'));    
    if(isempty(ppg_ind))
        ppg = [];
    else
        ppg = sig(:,ppg_ind);
        ppg(isnan(ppg)) = 0;
        signal=ppg;
    end    
      
    %CARGAMOS LOS RUIDOS
    
if(numRecord<=50)
    Noise = load('Noise_Seed_m03');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
   MA_Seed = [MA_Seed Rest];    
%    MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>50 && numRecord<=100)
    Noise = load('Noise_Seed_m06');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>100 && numRecord<=150)
    Noise = load('Noise_Seed_m08');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>150 && numRecord<=200)
    Noise = load('Noise_Seed_m12');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>200 && numRecord<=250)
    Noise = load('Noise_Seed_m13');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>250 && numRecord<=300)
    Noise = load('Noise_Seed_m15');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>300 && numRecord<=350)
    Noise = load('Noise_Seed_m17');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>350 && numRecord<=400)
    Noise = load('Noise_Seed_m1');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>400 && numRecord<=450)
    Noise = load('Noise_Seed_m053');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>450 && numRecord<=500)
    Noise = load('Noise_Seed_m07');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>500 && numRecord<=550)
    Noise = load('Noise_Seed_m09');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>550 && numRecord<=600)
    Noise = load('Noise_Seed_m11');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>600 && numRecord<=650)
    Noise = load('Noise_Seed_m115');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>650 && numRecord<=700)
    Noise = load('Noise_Seed_m14');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
if(numRecord>700)
    Noise = load('Noise_Seed_m085');
    MA_Seed = Noise.TotalGaussianNoise; 
    Rest = MA_Seed(1:(length(signal)-length(MA_Seed)));
    MA_Seed = [MA_Seed Rest];
%     MA_Seed = [MA_Seed zeros(1,length(signal(:,3))-length(MA_Seed))];
end
   
    %SUMAMOS A LA SEÑAL
    
    if(isempty(abp_ind))
        abp = [];
    else
        M=max(abp);
        m=min(abp);
        coef1=M-m;
        MA_Scaled=(MA_Seed*coef1)+m;
        abp=abp+MA_Scaled';
    end
    if(isempty(ppg_ind))
        ppg = [];
    else
        M=max(ppg);
        m=min(ppg);
        coef1=M-m;
        MA_Scaled=(MA_Seed*coef1)+m;
        ppg=ppg+MA_Scaled';
    end  

end

