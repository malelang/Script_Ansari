function [AUCNoisy, AUCOriginal, AUCDenoised,xCoordConvNoisy,xCoordConvOriginal,xCoordConvDenoised,yCoordConvNoisy,yCoordConvOriginal,yCoordConvDenoised]=generateConvexROC(reallabels,forecastOriginal,forecastNoisy,forecastDenoised,bradyflag)

[TPR_empOriginal, FPR_empOriginal,TPR_empNoisy,FPR_empNoisy,TPR_empDenoised, FPR_empDenoised]=compareEmpiricalcurve(reallabels,forecastOriginal,forecastNoisy,forecastDenoised,bradyflag);
[k1,v1]=convhull(FPR_empNoisy,TPR_empNoisy);
[k2,v2]=convhull(FPR_empOriginal,TPR_empOriginal);
[k3,v3]=convhull(FPR_empDenoised,TPR_empDenoised);
k1=k1(1:end-2);
k2=k2(1:end-1);
k3=k3(1:end-1);
xCoordConvNoisy=FPR_empNoisy(k1);
xCoordConvOriginal=FPR_empOriginal(k1);
xCoordConvDenoised=FPR_empDenoised(k1);
yCoordConvNoisy=TPR_empNoisy(k1);
yCoordConvOriginal=TPR_empOriginal(k1);
yCoordConvDenoised=TPR_empDenoised(k1);

AUCT=[v1 v2 v3]+0.5;
AUCNoisy=AUCT(1);
AUCOriginal=AUCT(2);
AUCDenoised=AUCT(3);
cols = [200 45 43; 37 64 180; 0 176 80; 0 0 0]/255;
if(bradyflag==1)
    figure
    plot(FPR_empNoisy(k1),TPR_empNoisy(k1),'-', 'color', cols(2,:),'linewidth',2),hold on, plot(FPR_empOriginal(k2),TPR_empOriginal(k2),'-', 'color', cols(1,:),'linewidth',2),hold on, plot(FPR_empDenoised(k3),TPR_empDenoised(k3),'-', 'color', cols(3,:),'linewidth',2),hold on, plot([0 1],[0 1],'k-'),...
        title('Bradycardia ROC curve fitting with Convexifying process'),xlabel('FPR'),ylabel('TPR'), legend(['Noisy results. AUC=' num2str(AUCT(1))],['Original results. AUC=' num2str(AUCT(2))],['Denoised results. AUC=' num2str(AUCT(3))])
else
    figure
    plot(FPR_empNoisy(k1),TPR_empNoisy(k1),'-', 'color', cols(2,:),'linewidth',2),hold on, plot(FPR_empOriginal(k2),TPR_empOriginal(k2),'-', 'color', cols(1,:),'linewidth',2),hold on, plot(FPR_empDenoised(k3),TPR_empDenoised(k3),'--', 'color', cols(3,:),'linewidth',2),hold on, plot([0 1],[0 1],'k-'),...
        title('Tachycardia ROC curve fitting with Convexifying process'),xlabel('FPR'),ylabel('TPR'), legend(['Noisy results. AUC=' num2str(AUCT(1))],['Original results. AUC=' num2str(AUCT(2))],['Denoised results. AUC=' num2str(AUCT(3))])
end

end