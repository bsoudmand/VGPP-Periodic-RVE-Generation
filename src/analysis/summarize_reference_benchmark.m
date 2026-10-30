function result = summarize_reference_benchmark(inputFile)
%SUMMARIZE_REFERENCE_BENCHMARK Recompute published baseline reachability table.
T=readtable(inputFile,'TextType','string');
required={'Method','TargetVF','ReachedTarget','FinalVF','ElapsedTime_s'};
assert(all(ismember(required,T.Properties.VariableNames)),'Missing benchmark columns.');
groups=findgroups(string(T.Method),T.TargetVF);
Method=splitapply(@(x) string(x(1)),string(T.Method),groups);
TargetAF=splitapply(@(x) x(1),T.TargetVF,groups);
Count=splitapply(@numel,T.TargetVF,groups);
SuccessRate=splitapply(@mean,double(T.ReachedTarget),groups);
MeanAchievedAF=splitapply(@mean,T.FinalVF,groups);
RuntimeMean_s=splitapply(@mean,T.ElapsedTime_s,groups);
RuntimeSD_s=splitapply(@std,T.ElapsedTime_s,groups);
result=table(Method,TargetAF,Count,SuccessRate,MeanAchievedAF,RuntimeMean_s,RuntimeSD_s);
end
