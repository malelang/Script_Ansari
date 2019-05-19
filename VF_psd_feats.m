function feats = VF_psd_feats(sig,Fs)

    sig = sig - medfilt1(sig,Fs);

    sig = sig - mean(sig);
    
    sig_st = abs(st(sig.*[linspace(0,1,50),ones(1,1900),linspace(1,0,50)]'));
    st_sum = sum(sig_st);
              
    medfs = [];
    wlen = 4*Fs;
    steplen = 1*Fs;
    for i = 1:steplen:length(sig)-wlen+1
        win = st_sum(i:i+wlen-1);
        
        pr = periodogram(win-mean(win));
        pr_cm = cumsum(pr)/sum(pr);

        medf = find(pr_cm(1:end-1)<=0.5 & pr_cm(2:end)>0.5);

        medfs = [medfs medf];
    end
    
    if(isempty(medfs))
        feats = nan;
    else
        feats = min(medfs);
    end

end