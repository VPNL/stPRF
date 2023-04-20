function df = st_compute(df,stim,params)

numFolds = df(1).numFolds;

for fold = 1 :numFolds
    fprintf(1,'[%s]:Compute TC Prediction (fold %d) \n',mfilename,fold);
    
    testSet = df(fold).test_set;
    fold_stim =  stim(:,:,testSet);
    
    x0     = df(fold).x0';
    y0     = df(fold).y0';
    sigma  = df(fold).sigma.major';
    exponent = df(fold).exponent';
    theta  = zeros(size(x0));
    
    if ~strcmp(df(fold).roi_name,'all') % if ROI is given
        if strcmp(params.analysis.viewType,'Inplane')
            coords = [1:length(x0)]';
        else
            coords = df(fold).coordsIndex;
            x0   =   x0(coords);
            y0    =  y0(coords);
            sigma =  sigma(coords);
            exponent = exponent(coords);
        end
    else % if whole brain
        coords = [1:length(x0)]';
    end
    
    
    %%
    n = numel(x0);
    if n == 1
        s = [1 2];
    elseif n < 1000
        s = [[1:ceil(n./10):n-2] n+1];
    else
        s = [[1:ceil(n./1000):n-2] n+1];
    end
    
    nChan = getChanNumber(params);
    prediction = zeros(size(fold_stim,2)/(params.analysis.temporal.tr*params.analysis.temporal.fs), ...
        n,nChan,size(fold_stim,3));
    drawnow;tic;
    for n=1:numel(s)-1
        if mod(n,ceil(numel(params.analysis.spatial.x0)/10)) == 0
            fprintf(1,'[%s]:(grid: %d/%d) \n',mfilename,n,numel(params.analysis.spatial.x0));
        end
        tmpParam = params;
        
        tmpParam.analysis.spatial.sigmaMajor      = sigma(s(n):s(n+1)-1);
        tmpParam.analysis.spatial.sigmaMinor      = sigma(s(n):s(n+1)-1);
        tmpParam.analysis.spatial.theta           = theta(s(n):s(n+1)-1);
        tmpParam.analysis.spatial.x0              = x0(s(n):s(n+1)-1);
        tmpParam.analysis.spatial.y0              = y0(s(n):s(n+1)-1);
        tmpParam.analysis.temporal.param.exponent = exponent(s(n):s(n+1)-1);
        
        % make predictions
        tmp = stPredictBOLDFromStim(tmpParam, fold_stim);
        
        % store
        prediction(:,s(n):s(n+1)-1,:,:) = gather(tmp.predBOLD);
        if ismember(n, round((1:10)/10* numel(s)-1)) % every 10% draw a dot
            fprintf(1,'.');drawnow;
        end
        clear tmp;
    end
    
    if size(fold_stim,3)> 1
         prediction = concatRuns(prediction);
    end
        prediction = permute(prediction, [2 1 3]);

%     if ismatrix(prediction)
%        tmp= zeros(1,size(prediction,1),size(prediction,2));
%        tmp(1,:,:) = prediction;
%        prediction = tmp;
%     end
    df(fold).test_pred = prediction;
    
    
    
    
end

end

