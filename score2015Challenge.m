%
%
% This script will calculate the statistics of your algorithm for each type
% of alarm, based on the true answer sheet. Your final score for the
% challenge will be a function of these statistics on the hidden test set.
%
%
% This script requires that you first run generateValidationSet.m
%
%
% Written by Ikaro Silva, January 26, 2015
%            Last Modified:
%

clear all;
fid=fopen('answers.txt','r');
if(fid ~= -1)
    ANSWERS=textscan(fid,'%s %d','Delimiter',',','EmptyValue',0); %Set empty values to FALSE alarms by default
    fclose(fid);
else
    error('Could not open users answer.txt for scoring. Run the generateValidationSet.m script and try again.')
end

fid=fopen(['ALARMS'],'r');
if(fid ~= -1)
    GOLD=textscan(fid,'%s %s %d','Delimiter',',');
    fclose(fid);
else
    error('Could not open challenge/ALARMS for scoring. Exiting...')
end

N=length(GOLD{1});
%Result columnes are: true positives, false positive, false negative, true
%negatives
RECORDS=GOLD{1};
ALARMS=GOLD{2};
ALARM_TYPES=unique(ALARMS);
NTYPES=length(ALARM_TYPES);
results=zeros(NTYPES,4);
GOLD_TRUTH=GOLD{3};
TPTach=[];
TNTach=[];
FPTach=[];
FNTach=[];
TPBrad=[];
TNBrad=[];
FPBrad=[];
FNBrad=[];
reallabelstachy=[];
reallabelsbrady=[];
forecastbrady=[];
forecasttachy=[];
%We do not assume that the Gold-standar and the Answers are sorted in the
%same order, so we search for the location of the individual records in
%ANSWER file
for n=1:N
    
    alarm_ind=strcmp(ALARMS{n},ALARM_TYPES);
    if(isempty(alarm_ind))
        error(['Unexpected alarm type: ' ALARMS{n} ' . Expected alarm types are: ' ALARM_TYPES{:} ])
    end
    
    rec_ind=strmatch(RECORDS{n},ANSWERS{1});
    if(isempty(rec_ind))
        warning(['Could not find answer for record: ' RECORDS{n} , ' setting it to a false alarm.'])
        this_answer=0;
    else
        this_answer=ANSWERS{2}(rec_ind);
    end
    if(this_answer ~=0)
        %Positive cases
        if(GOLD_TRUTH(n) == 1)
            %True positive
            results(alarm_ind,1)=results(alarm_ind,1)+1;
        else
            %False positive
            results(alarm_ind,2)=results(alarm_ind,2)+1 ;
        end
    else
        %Negative cases
        if(GOLD_TRUTH(n) == 1)
            %False negative
            results(alarm_ind,3)=results(alarm_ind,3)+1;
        else
            %True negative
            results(alarm_ind,4)=results(alarm_ind,4)+1;
        end
    end
    if(strcmp(ALARMS{n},'Tachycardia'))
        reallabelstachy=[reallabelstachy GOLD_TRUTH(n)];
        forecasttachy=[forecasttachy this_answer];
        if(this_answer ~=0)
        %Positive cases
            if(GOLD_TRUTH(n) == 1)
                %True positive
                TPTach=[TPTach 1];
                FPTach=[FPTach 0];
                TNTach=[TNTach 0];
                FNTach=[FNTach 0];
            else
                %False positive
                TPTach=[TPTach 0];
                FPTach=[FPTach 1];
                TNTach=[TNTach 0];
                FNTach=[FNTach 0];
            end
        else
            %Negative cases
            if(GOLD_TRUTH(n) == 1)
                %False negative
                TPTach=[TPTach 0];
                FPTach=[FPTach 0];
                TNTach=[TNTach 0];
                FNTach=[FNTach 1];
            else
                %True negative
                TPTach=[TPTach 0];
                FPTach=[FPTach 0];
                TNTach=[TNTach 1];
                FNTach=[FNTach 0];
            end
        end
    elseif(strcmp(ALARMS{n},'Bradycardia'))
        reallabelsbrady=[reallabelsbrady GOLD_TRUTH(n)];
        forecastbrady=[forecastbrady this_answer];
        if(this_answer ~=0)
        %Positive cases
            if(GOLD_TRUTH(n) == 1)
                %True positive
                TPBrad=[TPBrad 1];
                FPBrad=[FPBrad 0];
                TNBrad=[TNBrad 0];
                FNBrad=[FNBrad 0];
            else
                %False positive
                TPBrad=[TPBrad 0];
                FPBrad=[FPBrad 1];
                TNBrad=[TNBrad 0];
                FNBrad=[FNBrad 0];
            end
        else
            %Negative cases
            if(GOLD_TRUTH(n) == 1)
                %False negative
                TPBrad=[TPBrad 0];
                FPBrad=[FPBrad 0];
                TNBrad=[TNBrad 0];
                FNBrad=[FNBrad 1];
            else
                %True negative
                TPBrad=[TPBrad 0];
                FPBrad=[FPBrad 0];
                TNBrad=[TNBrad 1];
                FNBrad=[FNBrad 0];
            end
        end
    end
end

reallabelsbrady=double(reallabelsbrady);
reallabelstachy=double(reallabelstachy);
forecastbrady=double(forecastbrady);
forecasttachy=double(forecasttachy);
figure(1);
AUCBrady=GenerateROCCurve(reallabelsbrady,forecastbrady,1);
figure(2);
AUCTachy=GenerateROCCurve(reallabelstachy,forecasttachy,0);

total=sum(results,2);
nresults=results./repmat(total,[1 4]);
scores = sum(results(:,[1 4]),2)./(sum(results(:,[1 2 4]),2)+results(:,3)*5); 
gross=sum(results)/sum(sum(results));

for n=1:NTYPES
    indent=repmat(['\t'],[1 4-round(length(ALARM_TYPES{n})/8)]);
    str=[ALARM_TYPES{n} ':' indent 'TP: %1.3f\tFP: %1.3f \tFN: %1.3f\tTN: %1.3f\tScore: %1.3f\n'];
    fprintf(str,nresults(n,1),nresults(n,2),nresults(n,3),nresults(n,4),scores(n))
end

indent=repmat(['\t'],[1 4-round(length('Average')/8)]);
str=['Average:' indent 'TP: %1.3f\tFP: %1.3f \tFN: %1.3f\tTN: %1.3f\n'];
fprintf(str,mean(nresults(:,1)),mean(nresults(:,2)),mean(nresults(:,3)),mean(nresults(:,4)))

indent=repmat(['\t'],[1 4-round(length('Gross')/8)]);
str=['Gross:' indent 'TP: %1.3f\tFP: %1.3f \tFN: %1.3f\tTN: %1.3f\n'];
fprintf(str,gross(:,1),gross(:,2),gross(:,3),gross(:,4))
Sensibilities = 100*[nresults(1,1)/(nresults(1,1)+nresults(1,3)) nresults(2,1)/(nresults(2,1)+nresults(2,3)) nresults(3,1)/(nresults(3,1)+nresults(3,3)) nresults(4,1)/(nresults(4,1)+nresults(4,3)) nresults(5,1)/(nresults(5,1)+nresults(5,3))];
Especificities = 100*[nresults(1,4)/(nresults(1,2)+nresults(1,4)) nresults(2,4)/(nresults(2,2)+nresults(2,4)) nresults(3,4)/(nresults(3,2)+nresults(3,4)) nresults(4,4)/(nresults(4,2)+nresults(4,4)) nresults(5,4)/(nresults(5,2)+nresults(5,4))];
Precision= 100*[nresults(1,1)/(nresults(1,1)+nresults(1,2)) nresults(2,1)/(nresults(2,1)+nresults(2,2)) nresults(3,1)/(nresults(3,1)+nresults(3,2)) nresults(4,1)/(nresults(4,1)+nresults(4,2)) nresults(5,1)/(nresults(5,1)+nresults(5,2))];
Fmeasure= [(2*Precision(1)*Sensibilities(1))/(Precision(1)+Sensibilities(1)) (2*Precision(2)*Sensibilities(2))/(Precision(2)+Sensibilities(2)) (2*Precision(3)*Sensibilities(3))/(Precision(3)+Sensibilities(3)) (2*Precision(4)*Sensibilities(4))/(Precision(4)+Sensibilities(4)) (2*Precision(5)*Sensibilities(5))/(Precision(5)+Sensibilities(5))];

str2='Sensitivities:  %1.3f\t  %1.3f\t  %1.3f\t  %1.3f\t %1.3f\t \n';
str3='Specificities:  %1.3f\t  %1.3f\t  %1.3f\t  %1.3f\t %1.3f\t \n';
str4='Precision:  %1.3f\t  %1.3f\t  %1.3f\t  %1.3f\t %1.3f\t \n';
str5='Fmeasure:  %1.3f\t  %1.3f\t  %1.3f\t  %1.3f\t %1.3f\t \n';

fprintf(str2,Sensibilities(1),Sensibilities(2),Sensibilities(3),Sensibilities(4),Sensibilities(5));
fprintf(str3,Especificities(1),Especificities(2),Especificities(3),Especificities(4),Especificities(5));
fprintf(str4,Precision(1),Precision(2),Precision(3),Precision(4),Precision(5));
fprintf(str5,Fmeasure(1),Fmeasure(2),Fmeasure(3),Fmeasure(4),Fmeasure(5));