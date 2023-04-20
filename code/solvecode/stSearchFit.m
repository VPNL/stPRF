function view = stSearchFit(view, params, doDecimate, varargin)
% stSearchFit - find minimum for retinotopic model per voxel
%


%-----------------------------------
%--- input handling
%-----------------------------------
if ~exist('view','var') || isempty(view)
    error('Need view struct');
end;
if ~exist('params','var') || isempty(params)
    % See first if they are stored in the view struct
    params = viewGet(view,'rmParams');
    % if not loaded load them:
    if isempty(params)
        view = rmLoadParameters(view);
        params = viewGet(view,'rmParams');
    end
    % but allow ROI definitions to change
    if view.selectedROI == 0
        params.wData = 'all';
    else
        params.wData = 'roi';
    end
end
if ~exist('doDecimate','var') || isempty(doDecimate) || doDecimate<2
    if optimget(params.analysis.fmins.options,'MaxIter')<=0
        stage = 'fFit';
    else
        stage = 'sFit';
    end
    doDecimate = false;
else
    stage = 'sFit-sm';
end

% additional input handling
if nargin > 3
    addArg = varargin;
    if numel(addArg) == 1
        addArg=addArg{1};
    end
else
    addArg = [];
end

% parse command line inputs:
desc = [];
for n=1:2:numel(addArg)
    data = addArg{n+1};
    fprintf(1,'[%s]:Resetting %s\n',mfilename,addArg{n});
    switch lower(addArg{n})
        case {'desc'}
            desc = data;
        otherwise
            error('Unknown additional input');
    end
end

% for backward compatibility
if ~isfield(params.analysis.fmins,'options')
    params.analysis.fmins.options = optimset('fmincon');
    params.analysis.fmins.options = optimset(params.analysis.fmins.options,...
        'Display',params.analysis.fmins.Display,...
        'MaxIter',params.analysis.fmins.MaxIter,...
        'TolX',params.analysis.fmins.TolX,...
        'TolFun',params.analysis.fmins.TolFun);
end



%-----------------------------------
%--- loading data
%-----------------------------------
% get rmFile. This is the model definition that will start as a
% starting point for our search.
% stimfile = getAllFiles('./Stimuli','images_and_params*',1);

rmFile = viewGet(view,'rmFile');
if isempty(rmFile)
    rmFile    = fullfile(dataDir(view),[params.matFileName{end} '-gFit.mat']);
end

% load the entire grid when the grid is there.
if  strcmp(params.wData,'roi')
    checkpoint = strsplit(params.matFileName{end},'-');
    checkpoint{end} ='all';
    params.matFileName2 = strjoin(checkpoint,'-');
    rmFile2    = fullfile(dataDir(view),[params.matFileName2 '-gFit.mat']);
    if isfile(rmFile2)
        rmFile    = rmFile2;
    else
        rmFile    = fullfile(dataDir(view),[params.matFileName{end} '-gFit.mat']);
    end
end

fprintf(1,'[%s]:loading previous %s fit file \n',  mfilename, rmFile);
gridResult = load(rmFile,'model','df');
model = gridResult.model;
df    = gridResult.df;
fprintf(1,'[%s]:previous %s  fit file  loaded \n',  mfilename,rmFile);

%% Set Params

params = st_cleanParams(params);
% % 
% get re-assign params and grab defualt values 
params = getSpatialParams(params,3);
params = getTemporalParams(params);
nChan  = getChanNumber(params);
params = getExtraParams(params);
stim   = getStim(params);



%%
prevSkip = 1;


%%
for fold = 1:size(df,2)
    

    tmpFoldFile = fullfile(dataDir(view), ...
        [params.matFileName{1} '-fold' num2str(fold) '.mat' ] );

    if params.analysis.cv ==1
        checkDir = fullfile('./checkpoint',[params.matFileName{end},'-fold', num2str(fold)]);
    else
        checkDir = fullfile('./checkpoint',params.matFileName{end});
    end
    params.analysis.checkDir = checkDir;

    
    if ~isfolder(checkDir)
        mkdir(checkDir)
    end


    if prevSkip == 1
        % skip if precomputed results are there
        
        if isfile(tmpFoldFile)
            
            prevFOLD = load(tmpFoldFile);
            
            % get params
            switch lower(params.analysis.temporalModel)
                case {'glm','1ch-glm','spatial'}
                    fields = {'searchFit'};
                case {'dn','1ch-dcts','dn-st'}
                    fields = {'tau1', 'weight', 'tau2', 'nn', 'delay', 'searchFit'};
                case {'2ch','2ch-exp-sig'}
                    fields = {'tau_s', 'tau_ae', 'Lp', 'Kp', 'Kn','weight','searchFit'};
                case {'3ch-stln','cst'}
                    %                 fields = fieldnames(params.analysis.temporal.param);
                    %                 fields{end} = 'searchFit';
                    fields = {'tau_s', 'tau_t','searchFit'};
            end
            
            for ff =1:length(fields)
                df(fold).(fields{ff}) = [];
                
            end
            
            df(fold) = prevFOLD.fold_df;
            clear prevFOLD;
            fprintf('skipping previous fold %d \n',fold);
            continue
            
        end
    end
    
    % compute varexp if it is not already computed
    if ~isfield(df,'varexp')
        % update previous Beta if needed
        df(fold).varexp =[];
        % apply beta and re-calculate betas
        df    = st_applyBeta(df,fold);
    end
    
    % update model by giving df params
    model = updateModel(df,model,fold,params);
    

    %%
    % roi check
    switch lower(params.wData)
        case {'roi'}
            
            % if no roi is selected: select one
            if view.selectedROI == 0
                switch lower(view.viewType)
                    
                    case 'inplane'
                        % for inplanes default to gray matter
                        filename = 'gray.mat';
                        try
                            view   = loadROI(view,filename);
                        catch ME
                            error('[%s]:Cannot load ROI (%s).',mfilename,filename);
                            rethrow(ME)
                        end
                        
                        
                    otherwise
                        % otherwise ask
                        filename = getROIfilename(view);
                        view     = loadROI(view,filename);
                        
                end
            end
            ROIcoords = view.ROIs(view.selectedROI).coords;
            
        otherwise
            ROIcoords = [];
            % do nothing
    end
    
    
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
    
    % give some feedback so we know we are going
    vethresh = params.analysis.fmins.vethresh;
    if isempty(ROIcoords)
        fprintf(1,'[%s]:Processing voxels with variance explained >= %.2f\n',...
            mfilename,vethresh);
    else
        fprintf(1,'[%s]:Processing voxels with variance explained >= %.2f in ROI: %s\n',...
            mfilename,vethresh,view.ROIs(view.selectedROI).name);
    end
    drawnow;
    



    %%
    % go loop over slices
    for slice=loopSlices
        % put in number of data points. Right now this is the same as
        % size(data,1)
        
        %-----------------------------------
        % now we extract only the data from that slice and put it in a
        % temporary structure that will be modified throughout.
        %-----------------------------------
        s = rmSliceGet(model,slice);
        
        
        % The fitting uses fmincon which can only use type double (not single)
        % so convert model struct to double
        for n=1:numel(s)
            f=fieldnames(s{n});
            for n2=1:numel(f)
                if isnumeric(s{n}.(f{n2}))
                    s{n}.(f{n2}) = double(s{n}.(f{n2}));
                end
            end
        end
        
        
        %-----------------------------------
        % Find voxels (voxel>vethresh AND in ROI) that will be
        % processed.
        %-----------------------------------
        if isfield(model{1},'varexpfitprf')
            varexp = s{1}.varexpfitprf;
        else
            varexp   = 1-s{1}.rss./s{1}.rawrss;
        end
        if isempty(ROIcoords)
            wProcess = find(varexp>=vethresh);
            %check if there is searchfit already done, exclude
        else
            allcoords = viewGet(view,'coords', slice);
            if strcmpi('inplane', viewGet(view, 'viewType'))
                wProcess = find(varexp>=vethresh);

            else
                [tmp, wProcess] = intersectCols(allcoords,ROIcoords);
            end
            wProcess = wProcess(varexp(wProcess)>=vethresh);
        end
        
        if isempty(ROIcoords)
            gatherCheckpoints(view,params,wProcess);
        end

        searchFitComplete = zeros(size(model{1}.x0));
        searchFitComplete(:,wProcess) = 1;

        
        fprintf(1,'[%s]:Processing   %d voxels\n',...
            mfilename,length(wProcess));

  
        % if no voxel in the slice is valid, move on to the next slice
        if isempty(wProcess), continue; end
        
        
        
        %-----------------------------------
        % Place datasets behind each other. This is a rather crude way of
        % stimultaneously fitting both. Due to this we cannot
        % prewhiten (we could zeropadd/let the trends deal with this/not care).
        %-----------------------------------
        
        %%
        %-----------------------------------
        % Development Stage
        %-----------------------------------
        %     % Prepare data structure for cross-validation
        
        % we are doing cross validation here
        
        % perfom k-fold (3 fold in current example) cross valiation
%         rng('default') % For reproducibility
%         
         nStim = length(params.stim);
%         cv_split = cvpartition(nStim,'KFold',3);
%         numFolds = cv_split.NumTestSets;
        
        % Prepare data structure for cross-validation
        if params.analysis.cv == 1
            % perfom k-fold (3 fold in current example) cross valiation
            rng('default') % For reproducibility
            cv_split = cvpartition(nStim,'KFold',3);
            numFolds = cv_split.NumTestSets;
            
        else
            numFolds= 1;
            cv_split=[];
            
        end
     
        % create: data{fold}, prediction{fold}
        % train_grid (predictions) => Time X GRID X Channel
        % data => Time X Voxel
        % prediction and data should match in their time domain
   
        if params.analysis.cv
            trainSet =  find(cv_split.training(fold));
            testSet  =  find(cv_split.test(fold));
        else
            trainSet = 1:nStim;
            testSet = 1:nStim;
        end        % we get all the data
        
        p2 = params; 
        coarse   = false;
        if strcmpi( view.viewType, 'inplane' )
              p2.wData  = 'roi';
        else
            p2.wData = 'all';
        end
        
        [train_trend, train_ntrends, train_dcid] = rmMakeTrends(params,trainSet);
        [test_trend, test_ntrend, test_dcid] = rmMakeTrends(params,testSet);
        
        if ~isempty(params.analysis.userInputData) 
            if ~isfile(params.analysis.userInputData)
                error ('invalid input')
            end
            
            warning('loading **************userinputdata*************')
            load(params.analysis.userInputData)

            t1 = BOLD(:,:,trainSet);
            t1 = permute(t1,[1,3,2]);
            traindata_raw{fold} = reshape(t1,[],size(t1,3));
            
            t2 = BOLD(:,:,testSet);
            t2 = permute(t2,[1,3,2]);
            testdata_raw{fold} = reshape(t2,[],size(t2,3));
        else
            
            if params.analysis.coarseToFine ==1
                [traindata_raw{fold}, coords] = smoothData(view, p2, slice, ...
                    params.analysis.coarseToFine, [], trainSet);
                
                [testdata_raw{fold}, coords] = smoothData(view, p2, slice, ...
                    params.analysis.coarseToFine, [], testSet);
                
                % get ROI coords
                    if isempty(ROIcoords)
                        %wProcess = find(varexp>=vethresh);
                    else
                        coarseIndex = rmCoarseSamples(ROIcoords,params.analysis.coarseSample);
                    end
                
                % index into view's data

            else
                [traindata_raw{fold},~] = rmLoadData(view, p2, slice, ...
                    coarse, [], trainSet);
                [testdata_raw{fold},~] = rmLoadData(view, p2, slice, ...
                    coarse, [], testSet);
            end

        end
        
        if params.analysis.doDetrend
            
            trendBetas1 = pinv(train_trend)*traindata_raw{fold};
            trendBetas2 = pinv(test_trend)*testdata_raw{fold};
            
            traindata_raw{fold} = traindata_raw{fold} - train_trend*trendBetas1;
            testdata_raw{fold}  = testdata_raw{fold}  - test_trend*trendBetas2;

            trendBetas = trendBetas1;
        end
         
        % save cv information to a strcuct
        df(fold).info = cv_split;
        df(fold).numFolds = numFolds;
        df(fold).train_set = trainSet;
        df(fold).test_set = testSet;
        
        df(fold).train_data = [];
        df(fold).test_data =  [];
        
        df(fold).train_data_raw = traindata_raw{fold};
        df(fold).test_data_raw = testdata_raw{fold};
        
        df(fold).train_trend = train_trend;
        df(fold).test_trend = test_trend;
        
        df(fold).train_ntrends = train_ntrends;
        df(fold).test_ntrend = test_ntrend;
        
        df(fold).train_dcid = train_dcid;
        df(fold).test_dcid = test_dcid;
        

        %% SOLVE & Tidy and store data into 'df' structure
        % assign variables according to each fold
        data = df(fold).train_data_raw;
        trends = df(fold).train_trend; trends = single(trends);
        nTrends = df(fold).train_ntrends;
        dcid = df(fold).train_dcid;
        
       
        
        for n=1:numel(model)
            if isempty(desc)
                desc = lower(rmGet(model{n},'desc'));
            end
            switch desc
                case {'2d prf fit (x,y,sigma, positive only)','2',...
                        'difference 2d prf fit fixed (x,y,sigma,sigma2, center=positive)',...
                        'difference 2d prf fit beta fixed (x,y,sigma,sigma2, center=positive)',...
                        'oval 2d prf fit (x,y,sigma_major,sigma_minor,theta)', 'o',...
                        'radial oval 2d prf fit (x,y,sigma_major,sigma_minor)', 'r',...
                        'unsigned 2d prf fit (x,y,sigma)','u',...
                        'mirrored 2d prf fit (2*(x,y,sigma, positive only))','m',...
                        'shifted 2d prf fit (2*(x,y,sigma, positive only))',...
                        '1d prf fit (x,sigma, positive only)' ...
                        '2d nonlinear prf fit (x,y,sigma,exponent, positive only)', ...
                        '2d css nonlinear spatiotemporal prf fit', ...
                        '1ch spatiotemporal prf fit'
                        }
                    s{n}.b(2:nTrends+1,wProcess) = 0;
                    
                case {'double 2d prf fit (x,y,sigma,sigma2, center=positive)',...
                        'difference 2d prf fit (x,y,sigma,sigma2, center=positive)',...
                        'two independent 2d prf fit (2*(x,y,sigma, positive only))',...
                        'sequential 2d prf fit (2*(x,y,sigma, positive only))','s',...
                        'release two prf ties',...
                        'difference 1d prf fit (x,sigma, sigma2, center=positive)', ...
                        '2ch spatiotemporal prf fit'
                        }
                    s{n}.b(3:nTrends+2,wProcess) = 0;
                    
                otherwise
                    fprintf('Unknown pRF model: %s: IGNORED!',desc);
            end
        end
        
        
        % check to see that all t-series data contain finite numbers
        tmp      = sum(data(:,wProcess));
        ok       = ~isnan(tmp);
        wProcess = wProcess(ok); clear tmp ok;
        
        
        nSlices = length(loopSlices);

        for mm = 1:numel(model)
            model{mm} = rmSet(model{mm},'npoints',size(data,1));
            model{mm} = rmSet(model{mm},'pred_X', zeros(nSlices,size(data,2),model{mm}.npoints,nChan));
            for n=1:numel(s)
                s{n}.pred_X = zeros(model{mm}.npoints,size(data,2),nChan);
            end
        end
         
        
        
        
        % limit to voxels that will be processed
        data     = data(:,wProcess);

        % decimate? Note that decimated trends are stored in a new variable,
        data        = rmDecimate(data, doDecimate);
        sliceTrends = rmDecimate(trends, doDecimate);

        % store rawrss: this may be different from the one already there because
        % of the coarse-to-fine approach (i.e. smoothing). Please note that this
        % rawrss is the rss of the raw timeseries with the trends removed (i.e.
        % high-pass filtered.
        for n=1:numel(s)
            s{n}.rawrss(wProcess) = sum(double(data).^2);
        end
        
        % decimate predictions?
        %   If we have a nonlinear model, then we cannot pre-convolve the
        %   stimulus with the hRF. Instead we make predictions with the
        %   unconvolved images and then convolve with the hRF afterwards
        if ~checkfields(params, 'analysis', 'nonlinear') || ~params.analysis.nonlinear
            %�for a lineaer model, use the pre-convolved stimulus images
            original_allstimimages = params.analysis.allstimimages;
            params.analysis.allstimimages = rmDecimate(params.analysis.allstimimages,...
                doDecimate);
        else
            try
                % for a nonlinear model, use the unconvolved images
                params.analysis.allstimimages_unconvolved = rmDecimate(...
                    params.analysis.allstimimages_unconvolved, doDecimate);
                
                % scans stores the scan number for each time point. we need to keep
                % track of the scan number to ensure that hRF convolution does operate
                % across scans
                scans = rmDecimate(params.analysis.scan_number, doDecimate);
                params.analysis.scans = round(scans);
            end
        end
        
        %%
        %-----------------------------------
        % Go for each voxel
        %-----------------------------------
        
        for n=1:numel(model)
            % if dc is estimated from the data, remove it from the trends
            if params.analysis.dc.datadriven
                t.trends = [];
                t.dcid   = [];
            else
                t.trends = sliceTrends(:,dcid);
                t.dcid   = dcid;
            end
            if isempty(desc)
                desc = lower(rmGet(model{n},'desc'));
            end
            switch desc
                
                case {'st',  '1ch spatiotemporal prf fit', '2ch spatiotemporal prf fit'...
                        '2d css nonlinear spatiotemporal prf fit'}
                    
                    
                    fprintf('***************************************** \n');
                    fprintf('starting rmSearchFit_temporal_v3: \n');
                    fprintf('***************************************** \n');

                    trainSet = df(fold).train_set;
                    
                    % %% check gpu 
                    if params.useGPU
                        data = gpuArray(double(data));
                    else
                        data = double(data);
                    end
                        nChan = getChanNumber(params);
                        t.dcid           = t.dcid+nChan;

                        if isfield(params.analysis,'hrfver')
                            if params.analysis.hrfver > 0 
                                s{n}=rmSearchFit_temporal_v3_hrf(s{n},data,stim,params,wProcess,trainSet,t);
                            else 
                                s{n}=rmSearchFit_temporal_v3(s{n},data,stim,params,wProcess,trainSet,t);
                            end
                        else
                            s{n}=rmSearchFit_temporal_v3(s{n},data,stim,params,wProcess,trainSet,t);
                        end

                otherwise
                    fprintf('[%s]:Unknown pRF model: %s: IGNORED!\n',mfilename,desc);
            end
        end
        
        %%

    
        model = rmSliceSet(model,s,slice,nChan);

        
        % save and clean df output
        df(fold).x0            = model{1}.x0;
        df(fold).y0            = model{1}.y0;
        df(fold).sigma         = model{1}.sigma;
        df(fold).exponent      = model{1}.exponent;
        df(fold).beta          = model{1}.beta;
        df(fold).train_pred    = model{1}.pred_X;
        
        
        switch lower(params.analysis.temporalModel)
            case {'glm','1ch-glm','spatial'}
                
            case {'dn','1ch-dcts','dn-st'}
                df(fold).tau1          = model{1}.tau1;
                df(fold).tau2          = model{1}.tau2;
                df(fold).weight        = model{1}.weight;
                df(fold).nn            = model{1}.nn;
                df(fold).delay         = model{1}.delay;
            case {'2ch','2ch-exp-sig'}
                df(fold).tau_s         = model{1}.tau_s;
                df(fold).tau_ae        = model{1}.tau_ae;
                df(fold).Lp            = model{1}.Lp;
                df(fold).Kp            = model{1}.Kp;
                df(fold).Kn            = model{1}.Kn;
                df(fold).weight        = model{1}.weight;
                
            case {'3ch-stln','cst'}
                
                df(fold).tau_s        = model{1}.tau_s;
                df(fold).tau_t        = model{1}.tau_t;

        end

        
        
        % unlike grid-fit, varexp and train_pred is pre-computed
        df(fold).varexp        = model{1}.varexpfitprf;
               
        % remove a dimension
        px = model{1}.pred_X;
        px=reshape(px,[length(model{1}.x0),size(px,3),size(px,4)]);
        df(fold).train_pred    = sum((px),3)';
                
        % updateSearchFit index
        df(fold).searchFit = searchFitComplete;

        df(fold).roi_name      = 'all';
        df(fold).coords        = [];
        df(fold).coordsIndex   = [];
        if isfield(df,'train_grid'), df = rmfield( df , 'train_grid'); end
        if isfield(df,'test_grid'),  df = rmfield( df , 'test_grid' ); end

        fold_df = df(fold);
        save(tmpFoldFile,'fold_df','-v7.3')
    end    
end

% if there is no voxel
if ~isfield(df(1),'searchFit'), return; end
if params.analysis.hrfver == 3
    params.analysis.hrf.type = 'opt';
    [~,params] = getHRF(params);
    hrfs = params.analysis.hrf.func;    
end


%% perform crossValidation
%  - cv_varexp
%  - test_pred

% get stimulus if needed
if ~exist('stim','var')
    stim = getStim(params);
end

numFolds = df(1).numFolds;

for fold = 1:numFolds
    searchIDX = find(df(fold).searchFit);

    fprintf(1,'[%s]:Making %d model samples: ',mfilename,length(searchIDX));
for i =1:length(searchIDX)
    n = searchIDX(i);
    switch lower(params.analysis.temporalModel)
        case {'glm','1ch-glm','spatial'}
            fields ={};
        case {'dn','1ch-dcts','dn-st'}
            fields = {'tau1', 'weight','tau2','nn','delay'};
        case {'2ch','2ch-exp-sig'}
            error('need to do more work')
            %         fields = {'x0','y0','sigma','tau_s', 'tau_ae', 'Lp', 'Kp', 'Kn'};
        case {'3ch-stln','cst'}
            fields = {'exponent','tau_s', 'tau_t'};
    end

    
    tmpParam = params;
    tmpParam.analysis.spatial.x0         = df(fold).x0(n);
    tmpParam.analysis.spatial.y0         = df(fold).y0(n);
    tmpParam.analysis.spatial.sigmaMajor = df(fold).sigma.major(n);
    tmpParam.analysis.spatial.sigmaMinor = df(fold).sigma.minor(n);
    tmpParam.analysis.spatial.theta      = df(fold).sigma.theta(n);
    
    for nf = 1:length(fields)
        prm = df(fold).(fields{nf});
        if strcmp(fields{nf},'delay')
             fields{nf} = 'sigma';
        elseif  strcmp(fields{nf},'nn')
            fields{nf} = 'n';
        end
            tmpParam.analysis.temporal.param.(fields{nf}) = prm(n);
    end
        
    data =  df(fold).test_data_raw(:,n);      %    630       36559
    t.trends    =  df(fold).test_trend(:,df(fold).test_dcid);
    
    % set HRF
    if params.analysis.hrfver == 3
        tmpParam.analysis.hrf.func = hrfs(:,n);
    end

   
    % might want to use previous beta for re-fit
    [df(fold).test_pred(:,n), ~,  df(fold).cv_varexp(n)] ...
        = stReFit(tmpParam,stim(:,:,df(fold).test_set),data,t);

    if mod(i, 100) == 1 % every 10% draw a dot
        fprintf(1,'.');drawnow;
    end
end

    
end

%-----------------------------------
% save
%-----------------------------------
for n=1:numel(model)
    model{n} = rmSet(model{n},'coords',[]);
end
output = rmSave(view,model,params,1,stage,df);
view   = viewSet(view,'rmFile',output);

%-----------------------------------
% save and return output (if run interactively)
%-----------------------------------
% that's it
return;
%-----------------------------------


