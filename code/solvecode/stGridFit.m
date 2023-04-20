function view = stGridFit(view,params)
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



%%
[~, sDir ] =fileparts(fileparts(params.analysis.predFile));
if contains(sDir,'-')
    sDir2 = strsplit(sDir,'-');
    sDir2 = sDir2{1};
    params.analysis.predFile = strrep(params.analysis.predFile,sDir,sDir2);
end
if isfile(params.analysis.predFile)
    disp('***st predfile exists --- loading...')
    disp(params.analysis.predFile)
    load(params.analysis.predFile);  
elseif ~isfile(params.analysis.predFile)
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

% assign each grid to params
for es = 1:length(params.stim)
    params.stim(es).images_unconvolved= stim(:,:,es);
    params.stim(es).prediction = prediction(:,:,:,es);
end

      
        %        predBOLD: [210×35×2 gpuArray]
        %      predNeural: [21000×35×2 gpuArray]
        %          params: [1×1 struct]
        %            prfs: [3721×35 gpuArray]
        %     prfResponse: [21000×35×3 gpuArray]
        %             hrf: [2801×1 double]
        


%-----------------------------------
%--- now loop over slices
%--- but initiate stuff first
%-----------------------------------
switch lower(params.wData)
    case {'fig','roi'}
        loopSlices = 1;
    otherwise
        loopSlices = 1:params.analysis.nSlices;
end
nSlices = length(loopSlices);



for slice=loopSlices
    %-----------------------------------
    % Place datasets behind each other. This is a rather crude way of
    % stimultaneously fitting both. Due to this we cannot
    % prewhiten (we could zeropad/let the trends deal with this/not care).
    %-----------------------------------
    
    nStim = length(params.stim);
        
    % Prepare data structure for cross-validation
    if params.analysis.cv == 1
        % perfom k-fold (3 fold in current example) cross valiation
        rng('default') % For reproducibility
        cv_split = cvpartition(nStim,'KFold',3);
        numFolds = cv_split.NumTestSets;
    else
        %         prediction = [];
        numFolds= 1;        
        cv_split=[];
    end
    
    % create: data{fold}, prediction{fold}
    % train_grid (predictions) => Time X GRID X Channel
    % data => Time X Voxel
    % prediction and data should match in their time domain
    for fold = 1:numFolds
        if params.analysis.cv
            trainSet =  find(cv_split.training(fold));
            testSet  =  find(cv_split.test(fold));
        else
            trainSet = 1:nStim;
            testSet = 1:nStim;
        end
        % edited to account for cross validation
        [train_trend, train_ntrends, train_dcid] = rmMakeTrends(params,trainSet);
        [test_trend, test_ntrend, test_dcid] = rmMakeTrends(params,testSet);
        
        
        if ~isempty(params.analysis.userInputData) 
            if ~isfile(params.analysis.userInputData)
                error ('invalid input')
            end
            warning('**************userinputdata*************')
            load(params.analysis.userInputData)

            t1 = BOLD(:,:,trainSet);
            t1 = permute(t1,[1,3,2]);
            traindata{fold} = reshape(t1,[],size(t1,3));
            traindata_raw{fold} = traindata{fold};
            
            
            t2 = BOLD(:,:,testSet);
            t2 = permute(t2,[1,3,2]);
            testdata{fold} = reshape(t2,[],size(t2,3));
            testdata_raw{fold} = testdata{fold};
        else
            
            
            if params.analysis.coarseToFine
                [traindata{fold}, params,coords] = rmLoadData(view, params, slice, ...
                    params.analysis.coarseToFine, [], trainSet);
                
                [testdata{fold},params, coords] = rmLoadData(view, params, slice, ...
                    params.analysis.coarseToFine, [], testSet);
                
                [traindata_raw{fold}, coords] = smoothData(view, params, slice, ...
                    params.analysis.coarseToFine, [], trainSet);
                
                [testdata_raw{fold}, coords] = smoothData(view, params, slice, ...
                    params.analysis.coarseToFine, [], testSet);
                
                traindata{fold} = traindata_raw{fold};
                testdata{fold}  = testdata_raw{fold};

            else
                % concat and load  traning     data---
                [traindata{fold},params,~] = rmLoadData(view, params, slice, ...
                    params.analysis.coarseToFine, [], trainSet);
                [testdata{fold},params,~] = rmLoadData(view, params, slice, ...
                    params.analysis.coarseToFine, [], testSet);
                
                % concat and load  traning     data---
                [traindata_raw{fold},params,~] = rmLoadData(view, params, slice, ...
                    params.analysis.coarseToFine, [], trainSet);
                [testdata_raw{fold},params,~] = rmLoadData(view, params, slice, ...
                    params.analysis.coarseToFine, [], testSet);
            end
%             % save some raw unsmoothed values for future usage
%             [traindata_raw2{fold},params] = rmLoadData(view, params, slice, ...
%                 [], [], trainSet);
%             [testdata_raw{fold},params] = rmLoadData(view, params, slice, ...
%                 [], [], testSet);
        end
      
        
        % gather the prediction  data---
        train_grid = cell(length(trainSet),1);
        [train_grid{:}] = params.stim(trainSet).prediction;
%         train_grid = cell2mat(train_grid);
        
        test_grid = cell(length(testSet),1);
        [test_grid{:}] = params.stim(testSet).prediction;
        
        try
            train_grid = cell2mat(train_grid);
            test_grid = cell2mat(test_grid);
        catch 
            % do nothing
        end
        
        if params.analysis.doDetrend
            trendBetas1 = pinv(train_trend)*traindata{fold};
            trendBetas2 = pinv(test_trend)*testdata{fold};
            
            traindata{fold}     = traindata{fold} - train_trend*trendBetas1;
            testdata{fold}      = testdata{fold}  - test_trend*trendBetas2;
            
            trendBetas1 = pinv(train_trend)*traindata_raw{fold};
            trendBetas2 = pinv(test_trend)*testdata_raw{fold};
            
            traindata_raw{fold} = traindata_raw{fold} - train_trend*trendBetas1;
            testdata_raw{fold}  = testdata_raw{fold}  - test_trend*trendBetas2;
        end
        
        if params.analysis.doBlankBaseline
            
            tmp = cell(1,length(trainSet));
            [tmp{:}] = params.stim(trainSet).baseline;
            train_baselineIDX = cell2mat(tmp);
            
            tmp = cell(1,length(testSet));
            [tmp{:}] = params.stim(testSet).baseline;
            test_baselineIDX = cell2mat(tmp);
            clear tmp
            

            
            traindata2{fold}  = st_baselineCorrect(traindata{fold},train_baselineIDX);
            testdata{fold}  = st_baselineCorrect(testdata{fold},test_baselineIDX);
            
            traindata_raw{fold}  = st_baselineCorrect(traindata_raw{fold},train_baselineIDX);
            testdata_raw{fold}  = st_baselineCorrect(testdata_raw{fold},test_baselineIDX);
            
        end
        
        % save cv information to a strcuct
        df(fold).info = cv_split;
        df(fold).numFolds = numFolds;
        df(fold).train_set = trainSet;
        df(fold).test_set = testSet;
        
        df(fold).train_data = traindata{fold};
        df(fold).test_data = testdata{fold};
        
        df(fold).train_grid = train_grid;
        df(fold).test_grid = test_grid;
        
        df(fold).train_data_raw = traindata_raw{fold};
        df(fold).test_data_raw = testdata_raw{fold};
        
        df(fold).train_trend = train_trend;
        df(fold).test_trend = test_trend;
        
        df(fold).train_ntrends = train_ntrends;
        df(fold).test_ntrend = test_ntrend;
        
        df(fold).train_dcid = train_dcid;
        df(fold).test_dcid = test_dcid;
        
    end
    
    % remove variables that suck up storage resources
    clear train_grid;
    params.stim = rmfield( params.stim , 'prediction' ) ;
    params.stim = rmfield( params.stim , 'images_org' ) ;
    params.analysis = rmfield( params.analysis , 'allstimimages' ) ;
    params.analysis = rmfield( params.analysis , 'allstimimages_unconvolved' ) ;
end
%%
% SOLVE & Tidy and store data into 'df' structure
%     if params.analysis.cv == 1
numFolds = df(1).numFolds;
for fold = 1:numFolds
    
    % assign variables according to each fold
    
    data = df(fold).train_data;
    prediction = df(fold).train_grid;
    trends = df(fold).train_trend; trends = single(trends);
    ntrends = df(fold).train_ntrends;
    dcid = df(fold).train_dcid;
    
    % solve GRID! - for train_data set
    model = stGridSolve(params,data,prediction,trends,ntrends,dcid,slice,nSlices);
    
    % recreate complete model if we used coarse sampling
%     if params.analysis.coarseToFine
%         model = rmInterpolate(view, model, params);
%     end
    
    % setback to graymodel if it is an ROI based computation
    if isequal( lower(view.viewType),'gray' )
        model = rmGrayModel(view,model);
    end
    % save and clean df output
    df(fold).x0            = model{1}.x0;
    df(fold).y0            = model{1}.y0;
    df(fold).sigma         = model{1}.sigma;
    df(fold).exponent      = model{1}.exponent;
    df(fold).beta          = model{1}.beta;
    df(fold).train_pred    = model{1}.pred_X;
    
    
    
    switch lower(params.wData)
        case {'roi'}
            df(fold).roi_name      = params.roi.name;
            df(fold).coords        = params.roi.coords;
            df(fold).coordsIndex   = params.roi.coordsIndex;
        otherwise
            df(fold).roi_name      = 'all';
            df(fold).coords        = [];
            df(fold).coordsIndex   = [];
    end
%     clear model;
end

% save('ddf.mat','df','-v7.3')
% save("teest.mat","df","stim","params"d)
% apply beta and re-calculate betas
if ~isfield(df,'test_pred')
    %     df    = st_compute_test(df,params);
    df    = st_compute(df,stim,params);
end

% apply beta and calculate varexp
for fold = 1:numFolds
    df(fold).varexp =[];
    df    = st_applyBeta(df,fold,params);
    
    %             df(fold).roiVarExp = df(fold).varexp;
    %             df(fold).roiVarExp = df(fold).cv_varexp;
    
    if isequal( lower(view.viewType),'gray' )
%         model = rmGrayModel(view,model);
        model{1}.beta =[];
        model{1}.pred_X =[];
        model{1}.cv_varexpfitprf = df(fold).cv_varexp;
        model{1}.varexpfitprf = df(fold).varexp;
        Graymodel = rmGrayModel(view,model);
        df(fold).varexp = Graymodel{1}.varexpfitprf;
        df(fold).cv_varexp = Graymodel{1}.cv_varexpfitprf;
    end
    
    
end

% cleanup before saving
if isfield(df,'train_grid')
    df = rmfield( df , 'train_grid' );
end
if isfield(df,'test_grid')
    df = rmfield( df , 'test_grid' );
end



%-----------------------------------
% save and return output (if run interactively)
%-----------------------------------
rmFile = rmSave(view,model,params,1,'gFit',df);
view = viewSet(view,'rmFile',rmFile);
% saveSession;

% that's it
return;
%-----------------------------------


%-----------------------------------
%--- make all predictions first
%-----------------------------------
% nonlinear spatio-temproal model prediction-grid creation
% % if strcmp(params.analysis.pRFmodel{1}, 'st')
% %     
% %     % cache grid -- as it takes very very long time to generate
% %     % if we already created the gird, simply just load it!
% %     % stimGrid, prediction, grid are the same thing
% %     if exist(params.analysis.predFile, 'file')
% %         disp('***st predfile exists --- loading...')
% %         disp(params.analysis.predFile)
% %         load(params.analysis.predFile);
% % 
% %     elseif ~isfile(params.analysis.predFile)
% %         %         stimGrid = rmGridstPred(params);
% % %         stimGrid = rmGridstPred2(params);
% % %         stimGrid = rmGridstPred3(params);
% %         stimGrid = st_makeGrid(params);
% % 
% %     end
% %     % assign each grid to params
% %     for es = 1:length(params.stim)
% %         params.stim(es).prediction = single(full(stimGrid(es).prediction));
% %     end
% %     
% % %     prediction=cat(1,stimGrid.prediction);
% %     clear stimGrid;
% % else
% %     error("only works with spatiotemporal model")
% %     
% % end
% %%
% for es = 1:length(params.stim)
%     % check run-outputs and save
%     [tempPath,tempName]=fileparts(params.analysis.predFile);
%     runPred = strcat(tempName,'_r',num2str(es),'.mat');
%     cache_run_grid = fullfile(tempPath,runPred);
%     
%     if isfile(cache_run_grid)
%         load(cache_run_grid);read3dAnalyzeToTseries
%         stimGrid(es).prediction = prediction;
%         fprintf('Loading pre-saved Grid for %s model (stimulus = %d) \n', ...
%             params.analysis.temporalModel, es); drawnow;
%     elseif ~isfile(cache_run_grid)
%         n = numel(params.analysis.spatial.x0);
%         s = [[1:ceil(n./1000):n-2] n+1];
%         prediction = zeros(size(params.stim(es).images_unconvolved,2) ...
%             /(params.analysis.temporal.tr*params.analysis.temporal.fs), ...
%             n,length(unique(params.analysis.combineNeuralChan)));
%         fprintf(1,'[%s]:Making %d model samples:',mfilename,n);
%         drawnow;tic;
%         for n=1:numel(s)-1
%             if mod(n,ceil(numel(params.analysis.spatial.x0)/10)) == 0
%                 fprintf(1,'[%s]:(grid: %d/%d) \n',mfilename,n,numel(params.analysis.spatial.x0));
%             end
%             tmpParam = params;
%             
%             tmpParam.analysis.spatial.sigmaMajor = params.analysis.spatial.sigmaMajor(s(n):s(n+1)-1);
%             tmpParam.analysis.spatial.sigmaMinor = params.analysis.spatial.sigmaMinor(s(n):s(n+1)-1);
%             tmpParam.analysis.spatial.theta = params.analysis.spatial.theta(s(n):s(n+1)-1);
%             tmpParam.analysis.spatial.x0 = params.analysis.spatial.x0(s(n):s(n+1)-1);
%             tmpParam.analysis.spatial.y0 = params.analysis.spatial.y0(s(n):s(n+1)-1);
%             
%             % make rfs
%             tmp = stPredictBOLDFromStim(tmpParam, params.stim(es).images_unconvolved);
%             % store
%             prediction(:,s(n):s(n+1)-1,:) = tmp.predBOLD;
%             if ismember(n, round((1:10)/10* numel(s)-1)) % every 10% draw a dot
%                 fprintf(1,'.');drawnow;
%             end
%             clear tmp;
%         end
%         save(cache_run_grid, 'prediction','-v7.3')
%         stimGrid(es).prediction = prediction ;
%     end
% endf
% fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
% drawnow;

