%%
alarm_stats_creation_loader
% loaded 'ALARMS' file (converted it into csv format first
alarm_typs = unique(ALARMS.alarm_type);

% creating a table
alarm_stats = table(alarm_typs);

% total counts of each alarms
alarm_count_total = [];
for ii = 1:length(alarm_typs)
    alarm_count_total(ii,1) = numel(find(strcmp(ALARMS.alarm_type,alarm_typs(ii))));
end
alarm_stats = [alarm_stats,table(alarm_count_total)];
% count of truve/false alarms
true_alarms =[];
false_alarms =[];
for ii = 1:length(alarm_typs)
    actuals = ALARMS.actual(find(strcmp(ALARMS.alarm_type,alarm_typs(ii))));
    true_alarms(ii,1) = numel(find(actuals == 1));
    false_alarms(ii,1) = numel(find(actuals == 0));
end
alarm_stats = [alarm_stats,table(true_alarms)];
alarm_stats = [alarm_stats,table(false_alarms)];

%% checking quality of answers (loaded 'answers_entry' file)

% make it into a table
answers = table(ANSWERS{1,1},ALARMS.alarm_type,ALARMS.actual,ANSWERS{1,2},'VariableNames',{'file','alarm_typs','actuals','Output'}); 
alarm_typs = unique(answers.alarm_typs);

% creating a table
answers_stats = table(alarm_typs);

actuals = [];
predicted = [];
confmatrix = [];

for ii = 1:length(alarm_typs)
    fprintf('\n OUTPUT STATS ON %s \n',alarm_typs{ii});
    actuals(ii).act = answers.actuals(find(strcmp(ALARMS.alarm_type,alarm_typs(ii))));
    predicted(ii).pred = answers.Output(find(strcmp(ALARMS.alarm_type,alarm_typs(ii))));
    confmatrix(ii).conf = cfmatrix2(actuals(ii).act, predicted(ii).pred, [1 0], 0, 1);
 
% Outputs: confusion matrix
%
%                 Actual Classes
%                   p       n
%              ___|_____|______| 
%    Predicted  p'|     |      |
%      Classes  n'|     |      |

end
answers_stats = [answers_stats,table(confmatrix(:),actuals(:),predicted(:),'VariableNames',{'ConfusionMatrix','actuals','predicted'})];








