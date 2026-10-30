%% R1_C2_FixedDensityMC_Postprocess.m — MATLAB R2024b
% Reviewer 1 Comment 2: POSTPROCESS ONLY; never regenerates packings.
% Reads full VGPP particle coordinate CSVs saved by the efficient R1_C2 script.
% Hard-disk displacement Monte Carlo conserves N, r, periodic box and hard core.
% Samples along short chains show structural sensitivity to rearrangement;
% NOT a verified equilibrium hard-disk ensemble or causal proof of density-only order.
% Outputs CSV tables and two 300-dpi PNGs under revision_results/R1_C2_mc.
% All MC comparisons are paired to their own starting VGPP microstructure.
clc; clear; close all;

%% User settings
quickTestMode = false; % BEGIN TRUE. After smoke test, set FALSE.
fractions = [0.60 0.65 0.70];
seedBase = 9200; % Same as R1_C2_OrderingDecomposition_Efficient.m
nSeeds = 2;
nChainsPerSeed = 2; % independent MC chains initialized from same coordinate set
burnSweeps = 80; % N attempted moves/sweep
productionSweeps = 120;
nSnapshots = 6; % must divide productionSweeps
stepScale = 0.30; % trial displacement half-width / particle radius
boxLength = 10;
radiusFinal = 0.18;
overlapTolerance = 1e-8;

if quickTestMode
    fractions = 0.60;
    nSeeds = 1;
    nChainsPerSeed = 1;
    burnSweeps = 10;
    productionSweeps = 20;
    nSnapshots = 2;
end
assert(mod(productionSweeps,nSnapshots)==0,'nSnapshots must divide productionSweeps.');

% Can run from same folder as the prior script; no need to copy into cache.
resultRoot = 'revision_results';
cacheDir = fullfile(resultRoot,'R1_C2_efficient','coordinate_cache');
otherCacheDir = fullfile(resultRoot,'coordinate_cache');
outDir = fullfile(resultRoot,'R1_C2_mc');
if ~exist(outDir,'dir'), mkdir(outDir); end
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

%% Inspect cache, require at least one real saved configuration
requested = strings(0,1);
for iv=1:numel(fractions)
    for is=1:nSeeds
        % Efficient experiment's index: [0.60 0.65 0.70].
        allPhi = [0.60 0.65 0.70];
        ix = find(abs(allPhi-fractions(iv))<1e-9,1);
        seed = seedBase + 100*ix + is;
        requested(end+1,1)=sprintf('R1_C2_phi%02d_seed%d_Full_VGPP.csv',round(100*fractions(iv)),seed); %#ok<SAGROW>
    end
end
fprintf('Coordinate search: %s OR %s\n',cacheDir,otherCacheDir);
for k=1:numel(requested)
    fprintf('%s : %d\n',requested(k),isfile(fullfile(cacheDir,requested(k))) || isfile(fullfile(otherCacheDir,requested(k))));
end

%% Main MC: checkpoint after each chain
measureRows=table(); sweepRows=table(); missing=strings(0,1);
for iv=1:numel(fractions)
    phi=fractions(iv);
    allPhi=[0.60 0.65 0.70];
    ix=find(abs(allPhi-phi)<1e-9,1);
    nTarget=round(phi*boxLength^2/(pi*radiusFinal^2));
    for is=1:nSeeds
        seed=seedBase+100*ix+is;
        filename=sprintf('R1_C2_phi%02d_seed%d_Full_VGPP.csv',round(phi*100),seed);
        file=fullfile(cacheDir,filename);
        if ~isfile(file), file=fullfile(otherCacheDir,filename); end
        if ~isfile(file) && is==1
            legacyName=sprintf('R1_C2_%d_Full_VGPP_positions.csv',round(phi*100));
            legacyPath=fullfile(resultRoot,'R1_C2_efficient',legacyName);
            if ~isfile(legacyPath), legacyPath=fullfile(resultRoot,legacyName); end
            if isfile(legacyPath), file=legacyPath; end
        end
        if ~isfile(file)
            missing(end+1,1)=filename; %#ok<SAGROW>
            fprintf(2,'MISSING coordinates: %s\n',filename);
            continue;
        end
        T=readtable(file);
        assert(all(ismember({'x','y','radius'},T.Properties.VariableNames)),'Required columns x,y,radius missing.');
        x0=T.x(:); y0=T.y(:); r0=T.radius(:);
        assert(numel(x0)==nTarget,'Particle count mismatched for %s',filename);
        assert(all(abs(r0-radiusFinal)<1e-10),'Final radii not correct in %s',filename);
        assert(all(isfinite(x0)) && all(isfinite(y0)),'Invalid positions.');
        initialOverlap=R1C2MC_overlap_max(x0,y0,r0,boxLength);
        assert(initialOverlap<=overlapTolerance+1e-12,'Input has overlapping particles beyond tolerance: %s',filename);
        orig=R1C2MC_order(x0,y0,boxLength,radiusFinal);
        for chain=1:nChainsPerSeed
            rng(600000+1000*seed+chain,'twister');
            x=x0;y=y0;xu=x0;yu=y0;
            attempts=0;accepted=0;
            snapshots=burnSweeps+(productionSweeps/nSnapshots)*(1:nSnapshots);
            totalSweeps=burnSweeps+productionSweeps;
            for sweep=1:totalSweeps
                order=randperm(nTarget);
                batchAccepted=0;
                for j=1:nTarget
                    i=order(j);
                    dxStep=(2*rand-1)*stepScale*radiusFinal;
                    dyStep=(2*rand-1)*stepScale*radiusFinal;
                    tx=mod(x(i)+dxStep,boxLength);
                    ty=mod(y(i)+dyStep,boxLength);
                    dx=x-tx;dy=y-ty;
                    dx=dx-boxLength*round(dx/boxLength);
                    dy=dy-boxLength*round(dy/boxLength);
                    rr2=dx.^2+dy.^2;rr2(i)=inf;
                    attempts=attempts+1;
                    if all(rr2 >= (2*radiusFinal-overlapTolerance)^2)
                        x(i)=tx;y(i)=ty;
                        xu(i)=xu(i)+dxStep;yu(i)=yu(i)+dyStep;
                        accepted=accepted+1;
                        batchAccepted=batchAccepted+1;
                    end
                end
                % Cost-effective mixing diagnostics every 5 sweeps.
                if mod(sweep,5)==0 || any(snapshots==sweep)
                    travel=hypot(xu-x0,yu-y0)/(2*radiusFinal);
                    sweepRows=[sweepRows;table(phi,seed,chain,sweep, ...
                        batchAccepted/nTarget,accepted/attempts,mean(travel), ...
                        mean(travel.^2),'VariableNames', ...
                        {'NominalPhi','Seed','Chain','Sweep','SweepAcceptance', ...
                         'CumulativeAcceptance','MeanTravelDiameters','MeanSqTravelDiameters2'})]; %#ok<AGROW>
                end
                if any(snapshots==sweep)
                    m=R1C2MC_order(x,y,boxLength,radiusFinal);
                    maxOverlap=R1C2MC_overlap_max(x,y,r0,boxLength);
                    assert(maxOverlap<=overlapTolerance+1e-12,'MC created an invalid overlap.');
                    % Snapshots within a chain are correlated and must not be treated
                    % as independent realizations for statistical inference.
                    measureRows=[measureRows;table(phi,seed,chain,sweep, ...
                        nTarget,orig.LocalPsi6,m.LocalPsi6, ...
                        m.LocalPsi6-orig.LocalPsi6,orig.GlobalPsi6,m.GlobalPsi6, ...
                        orig.ContactCoord,m.ContactCoord,accepted/attempts, ...
                        mean(hypot(xu-x0,yu-y0))/(2*radiusFinal), maxOverlap, ...
                        'VariableNames',{'NominalPhi','Seed','Chain','Sweep','N', ...
                        'StartLocalPsi6','MCLocalPsi6','DeltaLocalPsi6', ...
                        'StartGlobalPsi6','MCGlobalPsi6','StartContactCoord', ...
                        'MCContactCoord','Acceptance','MeanTravelDiameters','MaxOverlap'})]; %#ok<AGROW>
                end
            end
            fprintf('phi %.2f seed %d chain %d: psi6 %.4f -> %.4f, acceptance %.3f\n', ...
                phi,seed,chain,orig.LocalPsi6,m.LocalPsi6,accepted/attempts);
            writetable(measureRows,fullfile(outDir,'R1_C2_MC_snapshots.csv'));
            writetable(sweepRows,fullfile(outDir,'R1_C2_MC_mixing.csv'));
        end
    end
end
if ~isempty(missing)
    writetable(table(missing,'VariableNames',{'MissingFilename'}),fullfile(outDir,'R1_C2_MC_missing_coordinates.csv'));
end
if isempty(measureRows)
    error('No configurations found. Copy original cached Full_VGPP coordinate CSVs into %s. No packing simulations were performed.',cacheDir);
end

%% Report initial vs final per independent seed (chains share starting coordinates)
last=measureRows(ismember(measureRows.Sweep,snapshots(end)),:);
summary=table();
for phi=unique(last.NominalPhi)'
    sel=last(abs(last.NominalPhi-phi)<1e-9,:);
    summary=[summary;table(phi,height(sel),numel(unique(sel.Seed)), ...
        mean(sel.StartLocalPsi6),mean(sel.MCLocalPsi6), ...
        mean(sel.DeltaLocalPsi6),std(sel.DeltaLocalPsi6), ...
        mean(sel.Acceptance),mean(sel.MeanTravelDiameters), ...
        'VariableNames',{'NominalPhi','NChains','NUniqueStartingPackings', ...
        'StartLocalPsi6Mean','LastMCLocalPsi6Mean','DeltaPsi6Mean', ...
        'DeltaPsi6SD_acrossChains_notIndependent','MeanAcceptance','MeanTravelDiameters'})]; %#ok<AGROW>
end
writetable(summary,fullfile(outDir,'R1_C2_MC_summary.csv'));

fig=figure('Color','w','Position',[100 100 1100 620]);hold on;box on;
for phi=unique(measureRows.NominalPhi)'
    sub=measureRows(abs(measureRows.NominalPhi-phi)<1e-9,:);
    sweeps=unique(sub.Sweep);
    v=arrayfun(@(s)mean(sub.DeltaLocalPsi6(sub.Sweep==s)),sweeps);
    plot(sweeps,v,'-o','LineWidth',2,'DisplayName',sprintf('Area fraction %.2f',phi));
end
yline(0,'k--','HandleVisibility','off');
xlabel('Monte Carlo sweep (N attempted moves)');ylabel('Change in mean local |\psi_6| from starting configuration');
legend('Location','best');set(gca,'FontSize',17,'LineWidth',1.8,'FontName','Times New Roman');
exportgraphics(fig,fullfile(outDir,'R1_C2_MC_OrderingChange.png'),'Resolution',300);close(fig);
fig=figure('Color','w','Position',[100 100 1100 620]);hold on;box on;
for phi=unique(sweepRows.NominalPhi)'
    sub=sweepRows(abs(sweepRows.NominalPhi-phi)<1e-9,:);
    sweeps=unique(sub.Sweep);
    v=arrayfun(@(s)mean(sub.MeanTravelDiameters(sub.Sweep==s)),sweeps);
    plot(sweeps,v,'-o','LineWidth',2,'DisplayName',sprintf('Area fraction %.2f',phi));
end
xlabel('Monte Carlo sweep (N attempted moves)');ylabel('Mean cumulative particle travel / diameter');
legend('Location','best');set(gca,'FontSize',17,'LineWidth',1.8,'FontName','Times New Roman');
exportgraphics(fig,fullfile(outDir,'R1_C2_MC_MixingTravel.png'),'Resolution',300);close(fig);
fprintf('\nDone: post-processing only. No packings regenerated.\n');
fprintf('Independent starting configurations: %d; MC chains are NOT independent packings.\n',numel(unique(measureRows.Seed)));

%% Local helpers (self-contained)
function maxO=R1C2MC_overlap_max(x,y,r,L)
    n=numel(x);maxO=0;
    for i=1:n-1
        dx=x(i+1:end)-x(i);dy=y(i+1:end)-y(i);
        dx=dx-L*round(dx/L);dy=dy-L*round(dy/L);
        maxO=max(maxO,max([0; r(i)+r(i+1:end)-hypot(dx,dy)]));
    end
end
function M=R1C2MC_order(x,y,L,r)
    n=numel(x);[SX,SY]=meshgrid(-1:1,-1:1);
    xx=repmat(x,9,1)+kron(SX(:),ones(n,1))*L;
    yy=repmat(y,9,1)+kron(SY(:),ones(n,1))*L;
    dt=delaunayTriangulation(xx,yy);E=edges(dt);
    central=@(ind) ind>4*n & ind<=5*n;
    E=[E(central(E(:,1)),:);fliplr(E(central(E(:,2)),:))];
    ii=E(:,1)-4*n;jj=mod(E(:,2)-1,n)+1;
    keep=ii~=jj;pairs=unique([ii(keep),jj(keep)],'rows');
    ii=pairs(:,1);jj=pairs(:,2);
    dx=x(jj)-x(ii);dy=y(jj)-y(ii);
    dx=dx-L*round(dx/L);dy=dy-L*round(dy/L);
    z=exp(1i*6*atan2(dy,dx));
    nBonds=accumarray(ii,1,[n 1],@sum,0);
    assert(all(nBonds>=3),'Delaunay neighbor count invalid.');
    psi=(accumarray(ii,real(z),[n 1],@sum,0)+1i*accumarray(ii,imag(z),[n 1],@sum,0))./nBonds;
    M.LocalPsi6=mean(abs(psi));M.GlobalPsi6=abs(mean(psi));
    nc=0;
    for i=1:n-1
        dx=x(i+1:end)-x(i);dy=y(i+1:end)-y(i);
        dx=dx-L*round(dx/L);dy=dy-L*round(dy/L);
        nc=nc+nnz(hypot(dx,dy)<=2.2*r);
    end
    M.ContactCoord=2*nc/n;
end
