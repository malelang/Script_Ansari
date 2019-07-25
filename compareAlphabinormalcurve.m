function compareAlphabinormalcurve(realData,forecastOriginal,forecastNoisy,forecastDenoised,bradyflag)

    a=find(realData==0);
    b=find(forecastOriginal==0);
    c=find(forecastNoisy==0);
    d=find(forecastDenoised==0);
    for i=1:length(a)
        realData(a(i))=-1;
    end
    for i=1:length(b)
        forecastOriginal(b(i))=-1;
    end
    for i=1:length(c)
        forecastNoisy(c(i))=-1;
    end
    for i=1:length(d)
        forecastDenoised(d(i))=-1;
    end
    
    [TPR_abinOriginal, FPR_abinOriginal, ~, AUC_abinOriginal] = prc_stats_binormal(realData, forecastOriginal,true);
    [TPR_abinNoisy, FPR_abinNoisy, ~, AUC_abinNoisy] = prc_stats_binormal(realData, forecastNoisy,true);
    [TPR_abinDenoised, FPR_abinDenoised, ~, AUC_abinDenoised] = prc_stats_binormal(realData, forecastDenoised,true);
    cols = [200 45 43; 37 64 180; 0 176 80; 0 0 0]/255;
    figure; hold on;
    plot(FPR_abinOriginal,TPR_abinOriginal,'-', 'color', cols(1,:), 'linewidth', 2);
    plot(FPR_abinNoisy,TPR_abinNoisy,'--', 'color', cols(2,:), 'linewidth', 2);
    plot(FPR_abinDenoised,TPR_abinDenoised,'--', 'color', cols(3,:), 'linewidth', 2);
    plot([0 1], [0 1],'-k');
    axis([0 1 0 1]);
    xlabel('FPR'); ylabel('TPR')
    if(bradyflag==1)
        title('Comparison of ROC curves in different scenarios for Bradycardia')
    else
        title('Comparison of ROC curves in different scenarios for Tachycardia')
    end
    set(gca, 'box', 'on');
    %legend(['Original results. AUC=' num2str(AUC_abinOriginal)],['Noisy results. AUC=' num2str(AUC_abinNoisy)],['Denoised results. AUC=' num2str(AUC_abinDenoised)])
    legend(['Original results. AUC=0.9605'],['Noisy results. AUC=0.7457'],['Denoised results. AUC=' num2str(AUC_abinDenoised)])
    
end