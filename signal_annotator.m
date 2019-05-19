function signal_annotator(files)

    % If input is provided, only annotate those files; otherwise, annotate
    % all.
    alarms = readtable('ALARMS','ReadVariableNames',0,'ReadRowNames',1);    
    if(nargin == 1)
        alarms = alarms(files,:);
    else
        del = [];
        for i = 1:size(alarms,1)
            if(~strcmp(alarms{i,1},'Ventricular_Flutter_Fib') || alarms{i,2}==1) %change the type of alarm here
                del = [del i];
            end
        end
        alarms(del,:) = []; % uncomment if you wanna only annotate the
        % alarms that are specified above
        
        % alarms = alarms(randperm(size(alarms,1)),:);   % uncomment if you
        % wanna annotate is random order
    end    
    
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
   
    det_win = 0.1; % If the guide is on, search in a neighborhood of det_win seconds to find the peak
    
    flag = 1;
    
    fi = 1;
    while(true)
%         if(strcmp(files{fi},'v623l'))
%             flag = 0;            
%         end
%         if(flag)
%             fi = fi + 1;
%             continue;
%         end
        
        disp(files{fi});        
        
        [~,signal,Fs,siginfo] = rdmat(files{fi});       
        
        if(isempty(dir([files{fi},'.ann'])))
            pk_locs = [];
            pk_types = [];            
        else
%             fi
%             fprintf('skipped %s',files{fi});
%             fi = fi+1;
%             continue; %Added by Belle

            temp = dlmread([files{fi},'.ann']);
            pk_locs = temp(:,1)';
            pk_types = temp(:,2)';
        end

        description=squeeze(struct2cell(siginfo));
        description=description(4,:);

        % Resample signal to 125Hz
        Fs=Fs(1);
        if Fs~=125
            signal=resample(signal,125,Fs);
            Fs=125;
        end
        win_len = round(det_win*Fs);

        xn = Fs*5*60; % alarm position
        x0 = xn-Fs*16+1; % 16s before the alarm
        signal = signal(x0:xn,:);
        
        [ecg1_orig, ecg2_orig, abp, ppg, Fs, names] = preprocess(signal, description, Fs);               
        
%         ecg1_orig = ecg1_orig - medfilt1(ecg1_orig,round(Fs/2));
%         ecg2_orig = ecg2_orig - medfilt1(ecg2_orig,round(Fs/2));
%         abp = abp - medfilt1(abp,round(Fs/2));
%         ppg = ppg - medfilt1(ppg,round(Fs/2));
        
%         [bl,al] = butter(2,[0.1,4]*2/Fs);
%         abp = filter(bl,al,abp);
%         ppg = filter(bl,al,ppg);
        
        gd = 1; % gd = 1 (positive guide), 0 (no guide), -1 (negative guide)

        h = figure('units','normalized','outerposition',[0 0 1 1]);        
        h1 = subplot(4,1,1);                
        update_fig(h1, ecg1_orig, pk_locs, pk_types, names{1}, [files{fi}, ' - ' , truth{fi}, ' ', alarm_types{fi}, ' - Guide ', num2str(gd), ' - Annotating...']);
        
        h2 = subplot(4,1,2);
        update_fig(h2, ecg2_orig, [], [], names{2}, []);
        
        h3 = subplot(4,1,3);
        update_fig(h3, abp, [], [], names{3}, []);
        
        h4 = subplot(4,1,4);
        update_fig(h4, ppg, [], [], names{4}, []);
             
        while(true)
            [x,~,button] = ginput(1);            
            if(~isempty(button) && (button==71 || button==103)) % if the button 'g' or 'G' is pressed
                gd = mod(gd+2,3)-1;
                update_fig(h1, ecg1_orig, pk_locs, pk_types, names{1}, [files{fi}, ' - ' , truth{fi}, ' ', alarm_types{fi}, ' - Guide ', num2str(gd), ' - Annotating...']);
                continue;
            end
            if(isempty(x)) % if 'return' key is pressed, exit the while loop
                break;
            end
            if(gca~=h1) % if user clicks on axes other than the top one, don't do anything
                continue;
            end          
            
            x = round(x);
            
            if(button == 1 || button == 2) % annotate peak (type 1 and 2)
                if(gd==1) % find a positive peak in the neighborhood
                    win = ecg1_orig(max(x-win_len,1):min(x+win_len,end));
                    [~,loc] = findpeaks(win);
                    [~,idx] = max(win(loc));
                    loc = loc(idx)+max(x-win_len,1)-1;
                elseif(gd==-1) % find a negative peak in the neighborhood (for reversed peaks)
                    win = -ecg1_orig(max(x-win_len,1):min(x+win_len,end));
                    [~,loc] = findpeaks(win);
                    [~,idx] = max(win(loc));
                    loc = loc(idx)+max(x-win_len,1)-1;
                else % insert a peak exactly where it was clicked
                    loc = x;
                end
                
                if(~isempty(loc)) % add the peak
                    pk_locs = [pk_locs, loc];
                    pk_types = [pk_types, button];
                    [pk_locs,ia] = unique(pk_locs);
                    pk_types = pk_types(ia);
                    
                    update_fig(h1, ecg1_orig, pk_locs, pk_types, names{1}, [files{fi}, ' - ' , truth{fi}, ' ', alarm_types{fi}, ' - Guide ', num2str(gd), ' - Annotating...']);
                end
            elseif(button == 3) % if right click, remove the closest peak
                pk_dist = abs(pk_locs-x);
                [~,idx] = min(pk_dist);
                pk_locs(idx) = [];
                pk_types(idx) = [];
                update_fig(h1, ecg1_orig, pk_locs, pk_types, names{1}, [files{fi}, ' - ' , truth{fi}, ' ', alarm_types{fi}, ' - Guide ', num2str(gd), ' - Annotating...']);
            end            
        end
        update_fig(h1, ecg1_orig, pk_locs, pk_types, names{1}, [files{fi}, ' - ' , truth{fi}, ' ', alarm_types{fi}, ' - Done - Press (S)ave, (F)orward or (B)ack']);
        
        % wait for commonds
        ch = '';
        while(isempty(ch) || (ch ~= 'f' && ch ~= 'b'))
            waitforbuttonpress;
            ch = get(gcf,'CurrentCharacter');
            if(lower(ch) == 'f') % go to the next subject
                fi = min(fi+1,length(files));
            elseif(lower(ch) == 'b') % go to the previous subject
                fi = max(fi-1,1);
            elseif(lower(ch) == 's') % save the annotations
                dlmwrite([files{fi} '.ann'],[pk_locs',pk_types']); % write the annotation into a file
                update_fig(h1, ecg1_orig, pk_locs, pk_types, names{1}, [files{fi}, ' - ' , truth{fi}, ' ', alarm_types{fi}, ' - Saved - Press (S)ave, (F)orward or (B)ack']);
            end
        end
        close(h);
    end

end

function update_fig(h, ecg, pk_locs, pk_types, xlbl, ttl)    
    % inputs: axis handle, signal, peak locations, peak types, text on the
    % X axis, text above the plot
    hold(h,'off');
    plot(h,ecg);
    axis(h,'tight');
    hold(h,'all');
    plot(h,pk_locs(pk_types==1),ecg(pk_locs(pk_types==1)),'r+');
    plot(h,pk_locs(pk_types==2),ecg(pk_locs(pk_types==2)),'k+');
    title(h,ttl);
    grid(h,'on'); grid(h,'minor');
    xlabel(h,xlbl);
end

