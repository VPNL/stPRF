function params = st_cleanParams(params)

%-----------------------------------
%--- For speed we do our computations in single precision.
%--- But we output in double (for compatibility).
%-----------------------------------
params.analysis.x0               = single(params.analysis.x0);
params.analysis.y0               = single(params.analysis.y0);
params.analysis.sigmaMajor       = single(params.analysis.sigmaMajor);
params.analysis.sigmaMinor       = single(params.analysis.sigmaMinor);
params.analysis.theta            = single(params.analysis.theta);
params.analysis.exponent         = single(params.analysis.exponent);
params.analysis.X                = single(params.analysis.X);
params.analysis.Y                = single(params.analysis.Y);
% params.analysis.allstimimages    = single(params.analysis.allstimimages);
params.analysis.sigmaRatio       = single(params.analysis.sigmaRatio);
params.analysis.sigmaRatioInfVal = single(params.analysis.sigmaRatioInfVal);
params.analysis.sigmaRatioMaxVal = single(params.analysis.sigmaRatioMaxVal);


params.analysis.spatial.x0 = params.analysis.x0;
params.analysis.spatial.y0 = params.analysis.y0;
params.analysis.spatial.sigmaMajor = params.analysis.sigmaMajor;
params.analysis.spatial.sigmaMinor = params.analysis.sigmaMinor;
params.analysis.spatial.theta = params.analysis.theta;

% tmp = params.analysis;
% tmp = rmfield(tmp,'dc');
% tmp = rmfield(tmp,'Hrf');
% tmp = rmfield(tmp,'x0');
% tmp = rmfield(tmp,'y0');
% tmp = rmfield(tmp,'sigmaMajor');
% tmp = rmfield(tmp,'sigmaMinor');
% 
% tmp = rmfield(tmp,'theta');
% tmp = rmfield(tmp,'X');
% tmp = rmfield(tmp,'Y');
% tmp = rmfield(tmp,'sigmaMinor');
% 
% cleanParams.analysis = tmp;
end

