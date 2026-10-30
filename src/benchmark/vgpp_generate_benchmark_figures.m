function vgpp_generate_benchmark_figures(params, plotData)
%VGPP_GENERATE_BENCHMARK_FIGURES Write figures after benchmark validation.
targetList = plotData.targetList;
sweepSummary = plotData.sweepSummary;
VGPPParticleTables = plotData.VGPPParticleTables;
VGPPReachedMatrix = plotData.VGPPReachedMatrix;
VGPPLocalRSDMatrix = plotData.VGPPLocalRSDMatrix;
VGPPSeedMatrix = plotData.VGPPSeedMatrix;
methods = unique(sweepSummary.Method,'stable');
resultsFolder = params.resultsFolder;
picsFolder = params.figuresFolder;
    %% ---------------- Figure R1: reachability ----------------
    fig1 = figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        errorbar(sweepSummary.TargetVF(idx), sweepSummary.MeanFinalVF(idx), ...
            sweepSummary.StdFinalVF(idx), '-o', 'LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    plot(targetList,targetList,'k--','LineWidth',1.2,'DisplayName','Target');
    xlabel('Target area fraction');
    ylabel('Achieved area fraction');
    legend('Location','best'); grid off;

    subplot(1,2,2); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        plot(sweepSummary.TargetVF(idx), ...
            sweepSummary.SuccessRate_percent(idx),'-o','LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    xlabel('Target area fraction');
    ylabel('Success rate (%)');
    ylim([0 105]); legend('Location','best'); grid off;

    save_figure_safely(fig1,fullfile(picsFolder,'Fig3_Target_VF_Reachability.png'));

    %% ---------------- Figure R2: equal-density spatial quality ----------------
    % RSA spatial metrics are shown only where RSA actually reached the target.
    fig2 = figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        xPlot = sweepSummary.TargetVF(idx);
        yPlot = 100*sweepSummary.MeanLocalVFRSD(idx);
        ePlot = 100*sweepSummary.StdLocalVFRSD(idx);
        if methods(im) == "RSA-PBC"
            ok = sweepSummary.SuccessRate_percent(idx) >= 95;
            xPlot = xPlot(ok); yPlot = yPlot(ok); ePlot = ePlot(ok);
        end
        errorbar(xPlot,yPlot,ePlot,'-o','LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    xlabel('Target area fraction');
    ylabel('Local area-fraction RSD (%)');
    legend('Location','best'); grid off;

    subplot(1,2,2); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        xPlot = sweepSummary.TargetVF(idx);
        yPlot = 100*sweepSummary.MeanLowVFCellFraction(idx);
        ePlot = 100*sweepSummary.StdLowVFCellFraction(idx);
        if methods(im) == "RSA-PBC"
            ok = sweepSummary.SuccessRate_percent(idx) >= 95;
            xPlot = xPlot(ok); yPlot = yPlot(ok); ePlot = ePlot(ok);
        end
        errorbar(xPlot,yPlot,ePlot,'-o','LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    xlabel('Target area fraction');
    ylabel('Low-area-fraction cells (%)');
    legend('Location','best'); grid off;

    save_figure_safely(fig2,fullfile(picsFolder,'Fig4_Spatial_Heterogeneity.png'));

    %% ---------------- Figure R3: ordering and runtime ----------------
    fig3 = figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1); hold on; box on;
    idxVG = sweepSummary.Method == "VGPP";
    plot(sweepSummary.TargetVF(idxVG),sweepSummary.MeanNNDCV(idxVG), ...
        '-o','LineWidth',1.8,'DisplayName','VGPP: NND CV');
    yyaxis right
    plot(sweepSummary.TargetVF(idxVG),sweepSummary.MeanRDFFirstPeak(idxVG), ...
        '-s','LineWidth',1.8,'DisplayName','VGPP: RDF first peak');
    yyaxis left
    xlabel('Target area fraction');
    ylabel('Nearest-neighbor distance CV');
    yyaxis right
    ylabel('RDF first peak');
    grid off;

    subplot(1,2,2); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        errorbar(sweepSummary.TargetVF(idx),sweepSummary.MeanTime_s(idx), ...
            sweepSummary.StdTime_s(idx),'-o','LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    xlabel('Target area fraction');
    ylabel('Elapsed time (s)');
    legend('Location','best'); grid off;

    % Runtime cannot be negative. Clip the displayed lower bound at zero,
    % while retaining the full positive mean + SD range.
    allTimeUpper = sweepSummary.MeanTime_s + sweepSummary.StdTime_s;
    finiteTimeUpper = allTimeUpper(isfinite(allTimeUpper));
    if ~isempty(finiteTimeUpper)
        ylim([0, 1.05*max(finiteTimeUpper)]);
    else
        ylim([0, inf]);
    end

    save_figure_safely(fig3,fullfile(picsFolder,'Fig6_Ordering_Runtime.png'));

    %% ---------------- Figure R4: strict geometric audit ----------------
    fig4 = figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        plot(sweepSummary.TargetVF(idx),sweepSummary.MeanMinClearance(idx), ...
            '-o','LineWidth',1.8,'DisplayName',char(methods(im)));
    end
    yline(-params.overlapAuditTolerance,'k--', ...
        'HandleVisibility','off');
    xlabel('Target area fraction');
    ylabel('Mean minimum clearance');
    legend('Location','best'); grid off;

    subplot(1,2,2); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        plot(sweepSummary.TargetVF(idx),sweepSummary.MeanOverlappingPairs(idx), ...
            '-o','LineWidth',1.8,'DisplayName',char(methods(im)));
    end
    xlabel('Target area fraction');
    ylabel('Mean overlapping pairs beyond tolerance');
    legend('Location','best'); grid off;

    save_figure_safely(fig4,fullfile(picsFolder,'Fig5_Periodic_Overlap_Audit.png'));

    %% ---------------- Figure 7: representative fixed-radius configurations ----------------
    % Select one actual successful realization at VF = 0.50, 0.60, and 0.70.
    % The representative realization is the run whose local-VF RSD is closest
    % to the mean local-VF RSD of all successful runs at that target.
    representativeTargets = [0.50 0.60 0.70];
    representativeTables = cell(numel(representativeTargets),1);
    representativeSeeds = nan(numel(representativeTargets),1);
    representativeRSD = nan(numel(representativeTargets),1);
    representativeExactVF = nan(numel(representativeTargets),1);

    for kk = 1:numel(representativeTargets)
        [~,itRep] = min(abs(targetList-representativeTargets(kk)));

        validRuns = find(VGPPReachedMatrix(itRep,:) & ...
            isfinite(VGPPLocalRSDMatrix(itRep,:)));

        if isempty(validRuns)
            warning('No successful VGPP realization found for target VF %.2f. Fig. 7 panel skipped.', ...
                representativeTargets(kk));
            continue;
        end

        meanRSDtarget = mean(VGPPLocalRSDMatrix(itRep,validRuns),'omitnan');
        [~,loc] = min(abs(VGPPLocalRSDMatrix(itRep,validRuns)-meanRSDtarget));
        selectedRun = validRuns(loc);

        representativeTables{kk} = VGPPParticleTables{itRep,selectedRun};
        representativeSeeds(kk) = VGPPSeedMatrix(itRep,selectedRun);
        representativeRSD(kk) = VGPPLocalRSDMatrix(itRep,selectedRun);

        TRep = representativeTables{kk};
        representativeExactVF(kk) = sum(pi*TRep.radius.^2)/(params.Lx*params.Ly);

        % Save the exact particle coordinates used in Fig. 7.
        fileTag = sprintf('%02d',round(100*representativeTargets(kk)));
        writetable(TRep,fullfile(resultsFolder, ...
            ['Fig7_VGPP_VF' fileTag '_representative_coordinates.csv']));
    end

    % Save a compact record documenting which realizations were selected.
    RepresentativeTargetVF = representativeTargets(:);
    RepresentativeExactVF = representativeExactVF;
    RepresentativeSeed = representativeSeeds;
    RepresentativeLocalVFRSD = representativeRSD;

    representativeSummary = table(RepresentativeTargetVF,RepresentativeExactVF, ...
        RepresentativeSeed,RepresentativeLocalVFRSD);

    writetable(representativeSummary,fullfile(resultsFolder, ...
        'Fig7_representative_selection_summary.csv'));

    % Draw the three actual accepted configurations.
    fig7 = figure('Color','w','Position',[100 100 1350 430]);

    panelLetters = {'(a)','(b)','(c)'};

    for kk = 1:numel(representativeTargets)
        ax = subplot(1,3,kk);
        hold(ax,'on');
        box(ax,'on');

        TRep = representativeTables{kk};

        if isempty(TRep)
            axis(ax,'off');
            text(ax,0.5,0.5,'No successful realization', ...
                'Units','normalized','HorizontalAlignment','center');
            continue;
        end

        plot_periodic_particle_configuration(ax,TRep,params.Lx,params.Ly);

        title(ax,sprintf('$\phi_A = %.5f$',representativeExactVF(kk)), ...
            'Interpreter','latex','FontWeight','bold','FontSize',16);

        text(ax,0.025,0.965,panelLetters{kk}, ...
            'Units','normalized','FontWeight','bold', ...
            'VerticalAlignment','top','FontSize',15);

        % Keep particle count comfortably inside the panel.
        text(ax,0.50,0.055, ...
            sprintf('N = %d',height(TRep)), ...
            'Units','normalized','HorizontalAlignment','center', ...
            'FontWeight','bold','FontSize',11, ...
            'BackgroundColor','w','Margin',1);
    end

    save_figure_safely(fig7,fullfile(picsFolder, ...
        'Fig7_Representative_FixedRadius_Configurations.png'));

    fprintf('\nFixed-radius target-VF study completed.\n');
    fprintf('CSV files: %s\n', resultsFolder);
    fprintf('PNG figures: %s\n', picsFolder);
end
