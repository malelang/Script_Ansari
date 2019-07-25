function comparePerfcurve(realData,forecastOriginal,forecastNoisy,forecastDenoised,bradyflag)

    [XROCorig,YROCorig,~,AUCorig]=perfcurve(realData,forecastOriginal,1);
    [XROCnoi,YROCnoi,~,AUCnoi]=perfcurve(realData,forecastNoisy,1);
    [XROCden,YROCden,~,AUCden]=perfcurve(realData,forecastDenoised,1);
    
    cols = [200 45 43; 37 64 180; 0 176 80; 0 0 0]/255;
    figure; hold on;
    plot(XROCorig,YROCorig,'-', 'color', cols(1,:),'linewidth',2),hold on,
    plot(XROCnoi,YROCnoi,'-', 'color', cols(2,:),'linewidth',2),hold on,
    plot(XROCden,YROCden,'--', 'color', cols(3,:),'linewidth',2),hold on,
    plot([0 1], [0 1],'-k');
    axis([0 1 0 1]);
    xlabel('FPR'); ylabel('TPR')
    if(bradyflag==1)
        title('Comparison of ROC curves in different scenarios for Bradycardia')
    else
        title('Comparison of ROC curves in different scenarios for Tachycardia')
    end
    set(gca, 'box', 'on');
    legend(['Original results. AUC=' num2str(AUCorig)],['Noisy results. AUC=' num2str(AUCnoi)],['Denoised results. AUC=' num2str(AUCden)])
    
end