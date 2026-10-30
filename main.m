% Main entry point for the fixed-final-radius VGPP/RSA-PBC benchmark.
clc; clear; close all;
repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'config'));
addpath(fullfile(repoRoot,'src','core'));
addpath(fullfile(repoRoot,'src','benchmark'));
addpath(fullfile(repoRoot,'src','analysis'));
addpath(fullfile(repoRoot,'src','validation'));
params = vgpp_config(repoRoot);
params.sourceFingerprint = vgpp_source_fingerprint(repoRoot);
if ~isfolder(params.resultsFolder), mkdir(params.resultsFolder); end
if ~isfolder(params.figuresFolder), mkdir(params.figuresFolder); end

% A full run retains the original six fractions and 20 seeds per method.
% For a smoke test, set targetVFSweepRuns=2 in config/vgpp_config.m.
vgpp_self_test(params);
rng(10);
[summary, plotData] = run_target_vf_sweep_block(params);
disp(summary);
if params.targetVFSweepRuns == 20 && params.runRSAInTargetSweep
    if ~isfile(params.referenceFile)
        error('VGPP:MissingReference', ...
            'Full benchmark validation requires the reference dataset: %s',params.referenceFile);
    end
    report = verify_reference_results( ...
        fullfile(params.resultsFolder,'Table_R1_fixedRadius_raw.csv'), ...
        params.referenceFile, params.resultsFolder, params.validationTolerance);
    if ~all(report.Pass)
        error('VGPP:ReferenceMismatch', ...
            'Numerical validation failed. Check outputs/benchmark/reference_comparison.csv before plotting.');
    end
end
vgpp_generate_benchmark_figures(params,plotData);
if params.targetVFSweepRuns == 20 && params.runRSAInTargetSweep
    plot_runtime_fig4(fullfile(params.resultsFolder,'Table_R1_fixedRadius_raw.csv'), ...
        params.figuresFolder, params.resultsFolder);
end
