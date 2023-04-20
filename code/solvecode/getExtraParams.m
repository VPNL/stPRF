function params = getExtraParams(params)
params.recomputePredictionsFlag = true;
params.analysis.zeroPadPredNeuralFlag = true;
params.analysis.reluFlag = true;

params.analysis.spatial.sparsifyFlag = 0;
params.saveDataFlag=0;
params.verbose = 0;

% normalization falg
params.analysis.normNeuralChan = 0;
params.analysis.normAcrossRuns = 0;

% pre-defind HRF at ms resolution
params.analysis.hrf.type = 'vista';
params.analysis.temporal.param.fs = 100; % fix it to pass input
params.analysis.temporal.fs = 100; % fix it to pass input
% params.analysis.temporal.param.fs = 1000; % fix it to pass input
% params.analysis.temporal.fs = 1000; % fix it to pass input

[~,params] = getHRF(params);


%% for noise-less
if isfield(params.analysis,'userInputData') && ~isempty(params.analysis.userInputData)
    if  contains(params.analysis.userInputData, {'pure', 'noise0', 'noiseless'}) || contains(pwd, {'pure', 'noise0', 'noiseless'})
        params.analysis.doDetrend = 0;
        params.analysis.doBlankBaseline = 0;
        params = getOptimParams(params,1);
        params.analysis.optim.nEval = 1;
    else
        params.analysis.doDetrend = 1;
        params.analysis.doBlankBaseline = 0;
        params = getOptimParams(params,1);
    end
    params.analysis.fmins.vethresh = -999; % compute all voxels
else
    % for general usage
    params.analysis.doDetrend = 1;
    params.analysis.doBlankBaseline = 0;
    params = getOptimParams(params,1);
end

params.analysis.optim.iter = 50;
params.analysis.optim.nEval = 3;

end