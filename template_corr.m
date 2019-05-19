for jj=1:length(ann_ppg)-1
% Calculate correlation coefficients based on the template length
    beatbegin=ann_ppg(jj);
    beatend=ann_ppg(jj+1);
    if beatend-beatbegin>3*Fs
        beatend=beatbegin+3*Fs;
    end
    period = ppg(beatbegin:beatend); % Selecting the waveform between onsets
    period_interp(jj,:) = interp1(linspace(1,500,length(period)),period,1:500);   % interpolating it to the same length signals
    
end

% creating mean template
templt = mean(period_interp,1);

% fiding corr coef with newly created averaged template
psqi=corr(period_interp',templt');
psqi(psqi(jj)<0) = 0;

% plotting mechanism
% figure;hold on;
% for jj = 1:size(period_interp,1)
%     plot(period_interp(jj,:),'b');
% end
% plot(templt,'r'); hold off;
