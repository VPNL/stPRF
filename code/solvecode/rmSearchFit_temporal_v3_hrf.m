function model = rmSearchFit_temporal_v3_hrf(model,data, stim, params, wProcess,trainSet,t)
% 2020/04 KIS: added spatiotemporal
% tic
warning('**************[ rmSearchFit_temporal_v3_hrf]wProcess*************')

% profile on;
% convert to double just in case
params.analysis.X = double(params.analysis.X);
params.analysis.Y = double(params.analysis.Y);


% check gpu
% if params.useGPU
%     stimtmp = gather(stim);
%     datatmp = gather(data);
%     g = gpuDevice();
%     reset(g);
%     stim = gpuArray(stimtmp(:,:,trainSet));
%     data = gpuArray(datatmp);
%     clear stimtmp datatmp;
% else
%     stim = stim(:,:,trainSet);
% end
stim = stim(:,:,trainSet);

%% set search range

model.s = model.s_major;
params.analysis.fmins.expandRange = 5;
expandRange = params.analysis.fmins.expandRange;
model.exp = model.exponent;
params.analysis.exp = params.analysis.exponent;

% 
[range, ~] = stSearchFit_range(params,model,data);

range.start = range.start([1:3,6],:);
range.lower = range.lower([1:3,6],:);
range.upper = range.upper([1:3,6],:);



switch params.analysis.temporalModel
    case '1ch-glm'
        temporal_searchRange =[];
    case '3ch-stLN'
        % 3 temporal params to solve:
        % 1) sustained delay 2) transient delay  3) exponent
          temporal_searchRange = [10 ;4 ; 100];
% 
    case '1ch-dcts'
        % 4 temporal params to solve:
        %  ["tau1", weight, "tau2", "n", "delay/sigma"]
           temporal_searchRange = [0.5  0   0.5    3   0.25; ...
                                   0.01  0  0.01   1  0.01; ...
                                   1     0     1   6   0.5];


end    



lower_p = []; upper_p =[]; init_p=[]; rstart=[]; rstep=[];
plb =[]; pub=[];
for ii = 1:numel(wProcess)
    vi = wProcess(ii);
    
    switch lower(params.analysis.temporalModel)
        case {'glm','1ch-glm','spatial'}
             init_p(ii,:) = [range.start(1:3,vi)'];
             lower_p(ii,:) = range.lower(1:3,vi)';
             upper_p(ii,:) = range.upper(1:3,vi)';
             plb(ii,:) =  lower_p(ii,:);
             pub(ii,:) =  upper_p(ii,:);

        case {'dn','1ch-dcts','dn-st'}
            init_p(ii,:) = [range.start(1:3,vi)' temporal_searchRange(1,:)];
            lower_p(ii,:) = [range.lower(1:3,vi)' temporal_searchRange(2,:)];
            upper_p(ii,:) = [range.upper(1:3,vi)' temporal_searchRange(3,:)];
            
            % fixed weight param to be zero
            % "tau1"    "weight"    "tau2"    "n"    "sigma"    
            plb(ii,:) =  lower_p(ii,:);
            pub(ii,:) =  upper_p(ii,:);
        case {'3ch','3ch-stln','cst'}
            init_p(ii,:) = [range.start(1:4,vi)' temporal_searchRange(1,:)];
            lower_p(ii,:) = [range.lower(1:4,vi)' temporal_searchRange(2,:)];
            upper_p(ii,:) = [range.upper(1:4,vi)' temporal_searchRange(3,:)];
            lower_p(ii,4)   = 0.1;
            plb(ii,:) =  lower_p(ii,:);
            pub(ii,:) =  upper_p(ii,:);

    end

    rstart(ii,:) = range.start(1:4,vi);
    rstep(ii,:)  = range.step(:,vi);
end

%% hrfs
if params.analysis.hrfver ==1
    load('./HRFindex_v1.mat');
    params.analysis.hrf.type = 'library1';
    [~,params] = getHRF(params);
    HRFindex = HRFindex(wProcess);
    hrflib = params.analysis.hrf.lib(:,1:2001);
    hrfs = zeros(2001,length(wProcess));
    for hh = 1:length(wProcess)
        hrfs(:,hh) = hrflib(HRFindex(hh),:);
    end
elseif params.analysis.hrfver ==2
    load('./HRFindex_v2.mat');
    params.analysis.hrf.type = 'library2';
    [~,params] = getHRF(params);
    HRFindex = HRFindex(wProcess);
    hrflib = params.analysis.hrf.lib(:,1:2001);
    hrfs = zeros(2001,length(wProcess));
    for hh = 1:length(wProcess)
        hrfs(:,hh) = hrflib(HRFindex(hh),:);
    end
elseif params.analysis.hrfver == 3
    params.analysis.hrf.type = 'opt';
    [~,params] = getHRF(params);
    hrfs = params.analysis.hrf.func;
    hrfs = hrfs(:,wProcess);
    
end


%% checkpoint
checkDir=params.analysis.checkDir;

if  params.analysis.hrfver ==1
    checkDir = [checkDir '_v1'];
elseif params.analysis.hrfver ==2
     checkDir = [checkDir '_v2'];
elseif params.analysis.hrfver ==3
    checkDir = [checkDir '_v3'];

end

mkdir(checkDir)
checkFile = fullfile(checkDir,'checkpoint.mat');

if isfile(checkFile)
    try
        backupParam = load(checkFile);
    catch
        backupParam = load(fullfile(checkDir,'checkpoint_backup.mat'));
    end

    copyfile(checkFile,fullfile(checkDir,'checkpoint_backup.mat'));

    estimatedParams = backupParam.estimatedParams;
    processIDX = find(isnan(estimatedParams(1,:)));
else
    estimatedParams = NaN(size(init_p))';
    processIDX =  1:numel(wProcess);
end

%%
opts.Display                    = 'off';         % Level of display ("iter", "notify", "final", or "off")';
opts.UncertaintyHandling = 0;
opts.MaxFunEvals = params.analysis.optim.iter;  

switch lower(params.analysis.temporalModel)
    case {'dn','1ch-dcts','3ch-stln'}
        nfirstEval = params.analysis.optim.nEval;
    otherwise
        nfirstEval = params.analysis.optim.nEval;
end

fprintf('[%d] Number of Voxels ... considering prevoius checkpoints \n',...
   numel(processIDX));

for ii = processIDX % 843
    params.analysis.hrf.func = hrfs(:,ii);
    opts.MaxFunEvals = params.analysis.optim.iter;  
    if mod(ii,100) == 0
        fprintf('[%d/%d] estimating:  --- %s ---\n',...
            ii,numel(wProcess),(params.analysis.temporal.model));
    end
    
    lb  =   lower_p(ii,:);          % Lower bounds
    ub  =   upper_p(ii,:);             % Upper bounds

    pl =   plb(ii,:) ;              % Plausible lower bounds
    pu =   pub(ii,:) ;                % Plausible upper bounds
    
    R2 =[]; tmpP = NaN([size(estimatedParams,1) nfirstEval]);
    for fi = 1:nfirstEval
        x0 = pl + (pu-pl).*rand(1,numel(pl));
        x0(1:3) = init_p(ii,1:3);

        try
            [tmpP(:,fi), ~,~, ~] = bads(@(x) ...
                solve_spatiotemporal(x,params, stim, data(:,ii),t), ...
                x0,lb,ub,pl,pu,[],opts);
        catch
            fprintf('lower bound initialization')
            x0(1:3) = lb(1:3);
            [tmpP(:,fi), ~,~, ~] = bads(@(x) ...
                solve_spatiotemporal(x,params, stim, data(:,ii),t), ...
                x0,lb,ub,pl,pu,[],opts);
        end
            
        % re-fit 
        params = setSeachParams(tmpP(:,fi),params);
        [~, ~, R2(fi)]=stReFit(params,stim,data(:,ii),t);

    end
    [~ ,mi]=max(R2);
    estimatedParams(:,ii) = tmpP(:,mi)';

    % save estimated params for every 50 voxels
    if mod(ii,50) == 0 || (processIDX(end) == ii)
        save(fullfile(checkDir,'checkpoint.mat'),'estimatedParams')
    end

    
end

%%
nChan = getChanNumber(params);
if params.useGPU
    pred = zeros(size(data),'gpuArray');
    beta = zeros(nChan,numel(wProcess),'gpuArray');
else
    pred = zeros(size(data));
    beta = zeros(nChan,numel(wProcess));
end
R2=zeros(1,numel(wProcess));

%% 
for ii = 1:numel(wProcess)
    params.analysis.hrf.func = hrfs(:,ii);
    params = setSeachParams(estimatedParams(:,ii),params);
    [  pred(:,ii), beta(:,ii),  R2(ii)]=stReFit(params,stim,data(:,ii),t);

end

%%
for ii = 1:numel(wProcess)

    vi = wProcess(ii);
    outParams = estimatedParams(:,ii);

    switch lower(params.analysis.temporalModel)
        case {'glm','1ch-glm','spatial'}
            model.x0(vi)         = outParams(1);
            model.y0(vi)         = outParams(2);
            model.s(vi)          = outParams(3);
            model.s_major(vi)    = outParams(3);
            model.s_minor(vi)    = outParams(3);
            model.s_theta(vi)    = 0;

            model.varexpfitprf(vi) = R2(ii);
            model.b(1,vi)      = gather(beta(:,ii));
            model.pred_X(:,vi,1)   = gather(pred(:,ii)) ;

        case {'dn','1ch-dcts','dn-st'}

            model.x0(vi)         = outParams(1);
            model.y0(vi)         = outParams(2);
            model.s(vi)          = outParams(3);
            model.s_major(vi)    = outParams(3);
            model.s_minor(vi)    = outParams(3);
            model.s_theta(vi)    = 0;

            model.tau1(vi)       = outParams(4);
            model.weight(vi)     = outParams(5);
            model.tau2(vi)       = outParams(6);
            model.nn(vi)          = outParams(7);
            model.delay(vi)      = outParams(8);

            model.varexpfitprf(vi) = R2(ii);
            model.b(1,vi)      = gather(beta(:,ii));
            model.pred_X(:,vi,1)   = gather(pred(:,ii)) ;

        case {'3ch','3ch-stln','cst'}
            %
            model.x0(vi)         = outParams(1);
            model.y0(vi)         = outParams(2);
            model.s(vi)          = outParams(3);
            model.s_major(vi)    = outParams(3);
            model.s_minor(vi)    = outParams(3);
            model.s_theta(vi)    = 0;
            
            model.exponent(vi)   = outParams(4);
            model.tau_s(vi)      = outParams(5);
            model.tau_t(vi)      = outParams(5);
            
            model.varexpfitprf(vi) = gather(R2(ii));
            model.b([1 2],vi)      = gather(beta(:,ii));
            model.pred_X(:,vi,1)   = gather(pred(:,ii)) ;
    end

end



% end time monitor
et  = toc;
if floor(et/3600)>0
    fprintf(1,'Done [%d hours].\n',ceil(et/3600));
else
    fprintf(1,'Done [%d minutes].\n',ceil(et/60));
end



end


