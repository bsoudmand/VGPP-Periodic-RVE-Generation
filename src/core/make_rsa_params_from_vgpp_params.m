function rsaParams = make_rsa_params_from_vgpp_params(params)
    rsaParams = struct();
    rsaParams.Lx = params.Lx;
    rsaParams.Ly = params.Ly;
    rsaParams.particleMode = params.particleMode;
    rsaParams.rMono = params.rMono;
    rsaParams.rMean = params.rMean;
    rsaParams.rStd = params.rStd;
    rsaParams.rMin = params.rMin;
    rsaParams.rMax = params.rMax;
    rsaParams.minGap = params.minGap;
    if isfield(params,'rsaMaxTrials')
        rsaParams.maxTrials = params.rsaMaxTrials;
    else
        rsaParams.maxTrials = 5000000;
    end
    if isfield(params,'activeTargetVF') && isfinite(params.activeTargetVF)
        rsaParams.targetVF = params.activeTargetVF;
    else
        rsaParams.targetVF = inf;
    end
    if isfield(params,'rsaMaxFailedTrials')
        rsaParams.maxFailedTrials = params.rsaMaxFailedTrials;
    else
        rsaParams.maxFailedTrials = 200000;
    end
end
