function [TPR_empOriginal, FPR_empOriginal,TPR_empNoisy, FPR_empNoisy,TPR_empDenoised, FPR_empDenoised]=compareEmpiricalcurve(realData,forecastOriginal,forecastNoisy,forecastDenoised,bradyflag)

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
    
    [TPR_empOriginal, FPR_empOriginal, ~, AUC_empOriginal] = prc_stats_empirical(realData, forecastOriginal);
    [TPR_empNoisy, FPR_empNoisy, ~, AUC_empNoisy] = prc_stats_empirical(realData, forecastNoisy);
    [TPR_empDenoised, FPR_empDenoised, ~, AUC_empDenoised] = prc_stats_empirical(realData, forecastDenoised);
    cols = [200 45 43; 37 64 180; 0 176 80; 0 0 0]/255;
    figure; hold on;
    plot(FPR_empOriginal,TPR_empOriginal,'-^', 'color', cols(1,:), 'linewidth', 2);
    plot(FPR_empNoisy,TPR_empNoisy,'--', 'color', cols(2,:), 'linewidth', 2);
    plot(FPR_empDenoised,TPR_empDenoised,'-', 'color', cols(3,:), 'linewidth', 2);
    plot([0 1], [0 1],'-k');
    axis([0 1 0 1]);
    xlabel('FPR'); ylabel('TPR')
    if(bradyflag==1)
        title('Comparison of ROC curves in different scenarios for Bradycardia')
    else
        title('Comparison of ROC curves in different scenarios for Tachycardia')
    end
    set(gca, 'box', 'on');
    legend(['Original results. AUC=' num2str(AUC_empOriginal)],['Noisy results. AUC=' num2str(AUC_empNoisy)],['Denoised results. AUC=' num2str(AUC_empDenoised)])
    
end