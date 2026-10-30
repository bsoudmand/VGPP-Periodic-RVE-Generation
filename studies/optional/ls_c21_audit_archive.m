function ls_c21_audit_archive(outputFolder,projectFolder,L,rFinal,tol)
% Read-only audit. Never changes the production raw, summary, coordinates or MATs.
assert(isfolder(outputFolder), ['LS archive folder not found: %s\n', ...
    'Set projectFolder or outputFolder in USER CONFIGURATION.'], outputFolder);
rawFile=fullfile(outputFolder,'R2_C21_optimized_raw.csv');
assert(isfile(rawFile),'Production raw CSV not found: %s',rawFile);
R=readtable(rawFile,'TextType','string');
required={'CaseID','Mode','TargetPhi','RealizablePhi','N','Status', ...
 'ExactTargetSuccess','AchievedPhi','FinalRadius','Runtime_s','Collisions'};
assert(all(ismember(required,R.Properties.VariableNames)), ...
       'Production raw CSV does not contain expected V5 schema.');
R.CaseID=string(R.CaseID);R.Mode=string(R.Mode);R.Status=string(R.Status);
if ~islogical(R.ExactTargetSuccess)
  R.ExactTargetSuccess=ismember(lower(string(R.ExactTargetSuccess)),["true","1"]);
end
% Audit production only; smoke and pilot are excluded from final benchmark.
P=R(R.Mode=="production",:);
assert(~isempty(P),'No production records found.');
% Refuse duplicate cases; the raw file should contain one latest state per case.
assert(numel(unique(P.CaseID))==height(P),'Duplicate production CaseID rows detected.');
runFolder=fullfile(projectFolder,'revision_results','R2_C21_consolidated_audit');
if ~isfolder(runFolder),mkdir(runFolder);end
stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
reportFolder=fullfile(runFolder,['audit_',stamp]);
if isfolder(reportFolder),reportFolder=[reportFolder,'_repeat'];end
mkdir(reportFolder);
fid=fopen(fullfile(reportFolder,'audit_manifest.txt'),'w');
assert(fid>0,'Could not create audit manifest.');
fprintf(fid,'R2-C21 consolidated archival audit\nMATLAB: %s\nArchive: %s\n',version,outputFolder);
fprintf(fid,'Raw CSV: %s\nGeometric tolerance: %.17g\n',rawFile,tol);
fprintf(fid,'Audits do not rerun dynamics or alter source archives.\n');
fprintf(fid,'Censored trajectories are not geometric failures.\n');
fclose(fid);
Nrec=height(P);
Geom_CaseID=P.CaseID;Geom_Status=P.Status;
Geom_File=strings(Nrec,1);Geom_Particles=nan(Nrec,1);
Geom_AreaFraction=nan(Nrec,1);Geom_MinGap=nan(Nrec,1);
Geom_OverlapPairs=nan(Nrec,1);Geom_Radius=nan(Nrec,1);
Geom_Pass=false(Nrec,1);Geom_AuditNote=strings(Nrec,1);
for k=1:Nrec
    name=char(P.CaseID(k));
    csv=fullfile(outputFolder,[name,'_coords.csv']);
    if ~isfile(csv)
        Geom_AuditNote(k)="Coordinate file missing; record not independently validated";
        continue;
    end
    Geom_File(k)=string(csv);
    C=readtable(csv,'VariableNamingRule','preserve');
    columns=lower(string(C.Properties.VariableNames));
    ix=find(columns=="x",1);iy=find(columns=="y",1);
    ir=find(ismember(columns,["radius","r"]),1);
    if isempty(ix)||isempty(iy)||isempty(ir)
        Geom_AuditNote(k)="Missing x/y/radius columns";continue;
    end
    x=double(C{:,ix});y=double(C{:,iy});rs=double(C{:,ir});
    n=numel(x);Geom_Particles(k)=n;
    if n~=P.N(k)||any(~isfinite([x;y;rs]))||any(rs<=0)|| ...
       any(x<0|x>=L|y<0|y>=L)
        Geom_AuditNote(k)="Invalid counts, coordinates or radii";continue;
    end
    Geom_Radius(k)=mean(rs);
    Geom_AreaFraction(k)=sum(pi*rs.^2)/L^2;
    worst=inf;bad=0;
    for i=1:n-1
        dx=abs(x(i+1:end)-x(i));dx=min(dx,L-dx);
        dy=abs(y(i+1:end)-y(i));dy=min(dy,L-dy);
        gap=hypot(dx,dy)-(rs(i)+rs(i+1:end));
        worst=min(worst,min(gap));bad=bad+sum(gap < -tol);
    end
    Geom_MinGap(k)=worst;Geom_OverlapPairs(k)=bad;
    sameRadii=all(abs(rs-rs(1))<1e-10);
    finalOK=~P.ExactTargetSuccess(k) || ...
        (abs(rs(1)-rFinal)<=1e-10 && ...
         abs(Geom_AreaFraction(k)-P.RealizablePhi(k))<1e-8);
    Geom_Pass(k)=bad==0 && sameRadii && finalOK && ...
        abs(Geom_AreaFraction(k)-P.AchievedPhi(k))<1e-8;
    if Geom_Pass(k),Geom_AuditNote(k)="PASS";
    else,Geom_AuditNote(k)="Mismatch/overlap; inspect saved coordinates";end
end
G=table(Geom_CaseID,Geom_Status,Geom_File,Geom_Particles, ...
 Geom_AreaFraction,Geom_Radius,Geom_MinGap,Geom_OverlapPairs, ...
 Geom_Pass,Geom_AuditNote, ...
 'VariableNames',{'CaseID','Status','CoordinateFile','Particles', ...
 'RecomputedPhi','MeanRadius','MinPeriodicGap','PairsBeyondTolerance', ...
 'Pass','AuditNote'});
writetable(G,fullfile(reportFolder,'R2_C21_coordinate_audit.csv'));
phiList=unique(P.TargetPhi);
Summary_TargetPhi=zeros(0,1);Summary_RealizablePhi=zeros(0,1);
Summary_Attempted=zeros(0,1);Summary_Completed=zeros(0,1);
Summary_Censored=zeros(0,1);Summary_MeanCompletedTime_s=zeros(0,1);
Summary_CompletedAudited=zeros(0,1);
for q=1:numel(phiList)
  idx=abs(P.TargetPhi-phiList(q))<1e-10;
  tmp=P(idx,:);idxRows=find(idx);
  done=tmp.ExactTargetSuccess & tmp.Status=="REACHED_RADIUS";
  censored=startsWith(tmp.Status,"CENSORED");
  times=tmp.Runtime_s(done);
  avg=NaN;if ~isempty(times),avg=mean(times);end
  Summary_TargetPhi(end+1,1)=phiList(q);
  Summary_RealizablePhi(end+1,1)=tmp.RealizablePhi(1);
  Summary_Attempted(end+1,1)=height(tmp);
  Summary_Completed(end+1,1)=sum(done);
  Summary_Censored(end+1,1)=sum(censored);
  Summary_MeanCompletedTime_s(end+1,1)=avg;
  Summary_CompletedAudited(end+1,1)=sum(done & Geom_Pass(idxRows));
end
S=table(Summary_TargetPhi,Summary_RealizablePhi,Summary_Attempted, ...
 Summary_Completed,Summary_Censored,Summary_MeanCompletedTime_s, ...
 Summary_CompletedAudited, ...
 'VariableNames',{'TargetPhi','RealizablePhi','Attempted','Completed', ...
 'ComputationallyCensored','MeanCompletedTime_s','CompletedWithPassingCoordinateAudit'});
writetable(S,fullfile(reportFolder,'R2_C21_recomputed_summary.csv'));
% Existing VGPP comparison is included only if a genuine archived table is found.
comparison=fullfile(projectFolder,'revision_results','R2_C21_VGPP_LS_final_comparison.csv');
if ~isfile(comparison)
  comparison=fullfile(projectFolder,'R2_C21_VGPP_LS_final_comparison.csv');
end
if isfile(comparison)
    copyfile(comparison,fullfile(reportFolder,'archived_VGPP_LS_comparison_unmodified.csv'));
end
fig=figure('Visible','off','Color','w','Position',[100 100 880 520]);
bar(categorical(string(S.TargetPhi)),[S.Completed,S.ComputationallyCensored],'stacked');
ylabel('Number of LS realizations');xlabel('Nominal particle area fraction');
legend({'Reached final radius','Computationally censored'},'Location','best');
title('LS production outcomes (archived runs)');
set(gca,'FontName','Times New Roman','FontSize',13);
exportgraphics(fig,fullfile(reportFolder,'R2_C21_completion_summary.png'),'Resolution',250);
close(fig);
report=fopen(fullfile(reportFolder,'R2_C21_interpretation.txt'),'w');
fprintf(report,'No dynamics were executed. All outputs were derived from archived results.\n');
fprintf(report,'Completed coordinate audits passed: %d/%d available completion records.\n', ...
    sum(Geom_Pass & P.ExactTargetSuccess),sum(P.ExactTargetSuccess));
fprintf(report,'A missing coordinate file is not proof of an invalid LS configuration.\n');
fprintf(report,'Only completed runtimes enter reported means; censored run costs are excluded.\n');
fprintf(report,'Unequal VGPP/LS realization counts and implementation details prevent general speed ranking.\n');
for q=1:height(S)
 fprintf(report,'phi=%.2f completed=%d attempted=%d censored=%d meanCompletedTime=%.6g s\n', ...
 S.TargetPhi(q),S.Completed(q),S.Attempted(q),S.ComputationallyCensored(q),S.MeanCompletedTime_s(q));
end
fclose(report);
disp(S);
fprintf('Read-only R2-C21 audit saved to: %s\n',reportFolder);
end
