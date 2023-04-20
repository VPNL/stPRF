function params=getOptimParams(params,doRidge)

params.analysis.optim.decimate = 1; %decimate by 1 so, no decimate
params.analysis.optim.fprec = 1e-4;

if doRidge == 1
    params.analysis.optim.ridge = 1;
    params.analysis.optim.ridgeAlpha = 0;
    params.analysis.optim.algorithm = 'interior-point';

elseif doRidge == 0
    params.analysis.optim.ridge = 0;
    params.analysis.optim.ridgeAlpha = 0;
    params.analysis.optim.algorithm = 'interior-point';

end
% params.analysis.optim.algorithm ='sqp'
end