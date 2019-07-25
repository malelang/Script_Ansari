function[vectorROC,metrics]=GenerateROCinfo(results,flagBrad)
    [a,b]=size(results);
    cont=0;
    Resultsorganized=[];
    for i=1:a
        if(results(i,1)==1)
            %true positive, must go in first
            Resultsorganized(cont+1,:)=results(i,:);
            cont=cont+1;
        end
    end
    for i=1:a
        if(results(i,2)==1 || results(i,4)==1)
            %false positive or false negative
            Resultsorganized(cont+1,:)=results(i,:);
            cont=cont+1;
        end
    end
    for i=1:a
        if(results(i,3)==1)
            %true negative
            Resultsorganized(cont+1,:)=results(i,:);
            cont=cont+1;
        end
    end
    
    %See how many TP the results have to start from there
    %Once the vectors for FP, FN, TN and TP are organized again, the info
    %is retrieved for the ROC chart
    %For bradycardia results, we'll start from 44 and then 10 each time
    
    if(flagBrad)
        register=44;
        cont2=1;
        vectorROC(cont2,:)=sum(Resultsorganized((cont2:register),:));
        cont2=cont2+1;
        for i=register:10:a-10
            vectorROC(cont2,:)=sum(Resultsorganized((1:i+10),:));
            cont2=cont2+1;
        end 
        vectorROC(cont2,:)=sum(Resultsorganized((1:end),:));
    else
        register=100;
        cont2=1;
        vectorROC(cont2,:)=sum(Resultsorganized((cont2:register),:));
        cont2=cont2+1;
        for i=register:2:a-2
            vectorROC(cont2,:)=sum(Resultsorganized((1:i+2),:));
            cont2=cont2+1;
        end 
    end
    
    tpr=vectorROC(:,1)./(vectorROC(:,1)+vectorROC(:,4));
    tnr=vectorROC(:,3)./(vectorROC(:,2)+vectorROC(:,3));
    fnr=1-tnr;
    metrics=[tpr tnr fnr];
    metrics(end+1,:)=[0 1 0];
    
end