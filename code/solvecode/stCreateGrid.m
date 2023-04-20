function view = stCreateGrid(view,params)
% stGridFit - fit spatiotemporal retinotopic model along grid (coarse stage)
%
% the structure of the code follows rmGridFit
% 2020/12 Insub Kim:

if notDefined('view'),   error('Need view struct'); end
if notDefined('params'), error('Need params'); end


%-----------------------------------
%--- For speed we do our computations in single precision.
%--- But we output in double (for compatibility).
%-----------------------------------

params = st_cleanParams(params);

% get re-assign params and grab defualt values
params = getSpatialParams(params,3);
params = getTemporalParams(params);
nChan  = getChanNumber(params);
params = getExtraParams(params);
stim   = getStim(params);



% %% downsample stimulus and reshape according to the temporal sampling rate
% downSampleRate = 1000/params.analysis.temporal.param.fs;
% stim_all = [];
% if params.analysis.temporal.param.fs ~= 1000
%     for es = 1:length(params.stim)
%         stim = full(logical( params.stim(es).images_unconvolved'));
%         stim_all(:,:,es) = single(downsample(stim,downSampleRate)'); %3721 21000 6
%     end
% end
% clear stim;
% %% check gpu
% if params.useGPU
%     g = gpuDevice();
%     reset(g);
%     if canUseGPU()
%         stim = gpuArray(stim_all);
%     end
% else
%    stim = stim_all;
% end
% clear stim_all;
%%

n = numel(params.analysis.spatial.x0);
s = [[1:ceil(n./1000):n-2] n+1];
prediction = zeros(size(stim,2)/(params.analysis.temporal.tr*params.analysis.temporal.fs), ...
    n,nChan,size(stim,3));
fprintf(1,'[%s]:Making %d model samples: ',mfilename,n);
drawnow;tic;
for n=1:numel(s)-1
    tmpParam = params;
    tmpParam.analysis.spatial.sigmaMajor = params.analysis.spatial.sigmaMajor(s(n):s(n+1)-1);
    tmpParam.analysis.spatial.sigmaMinor = params.analysis.spatial.sigmaMinor(s(n):s(n+1)-1);
    tmpParam.analysis.spatial.theta = params.analysis.spatial.theta(s(n):s(n+1)-1);
    tmpParam.analysis.spatial.x0 = params.analysis.spatial.x0(s(n):s(n+1)-1);
    tmpParam.analysis.spatial.y0 = params.analysis.spatial.y0(s(n):s(n+1)-1);
    tmpParam.analysis.temporal.param.exponent = params.analysis.exponent(s(n):s(n+1)-1);
    % make rfs
    tmp = stPredictBOLDFromStim(tmpParam, stim);
    % store
    prediction(:,s(n):s(n+1)-1,:,:) = tmp.predBOLD;
    if ismember(n, round((1:10)/10* numel(s)-1)) % every 10% draw a dot
        fprintf(1,'.');drawnow;
    end
    clear tmp;
end
save(params.analysis.predFile, 'prediction','-v7.3')
fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
drawnow;

end

