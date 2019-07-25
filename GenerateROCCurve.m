function AUC=GenerateROCCurve(reallabels,forecastlabels,flagbrady)
    
    [XROC,YROC,~,AUC]=perfcurve(reallabels,forecastlabels,1);
    plot(XROC,YROC),xlabel('False Positive Rate'),ylabel('True Positive Rate')
    if(flagbrady==1)
        title('ROC Curve for Bradycardia Classification with Sardar Ansari Script')
    else
        title('ROC Curve for Tachycardia Classification with Sardar Ansari Script')
    end
    
end