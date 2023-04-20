function df = st_applyBeta(df,fold,params)

% recompute and clear unwanted-variables
% 1) remove grid info (no longer needed)
% 2) compute & override trainPred with train_pred * beta
% 3) calculate varexp

if isfield(df,'test_pred')
    makeTestPred = 1;
else
    makeTestPred = 0;
end

if isempty(df(fold).varexp) || ndims(df(1).train_pred) == 4
    
    x0     = df(fold).x0';
    y0     = df(fold).y0';
    sigma  = df(fold).sigma.major';
    exponent = df(fold).exponent';
    trainData_raw = df(fold).train_data_raw;
    testData_raw = df(fold).test_data_raw;
    
    
    % squeeze only the first dim
    V = df(fold).train_pred;
    sz = size(V); t2 = sz~=1; t2(2)=1;
    pred_X_train = reshape(V,sz(t2)); 
    
    if makeTestPred == 1
        V = squeeze(df(fold).test_pred);
        if ndims(V) == 2
            pred_X_test = df(fold).test_pred;
        else
            pred_X_test = V;
        end
    end
    
    % if ROI is given
    if ~strcmp(df(fold).roi_name,'all') && ~strcmp(params.analysis.viewType,'Inplane')
        coords = df(fold).coordsIndex;
        x0   =   x0(coords);
        y0    =  y0(coords);
        sigma =  sigma(coords);
        exponent = exponent(coords);
        V = df(fold).beta(:,coords,:);
        sz = size(V); t2 = sz~=1; t2(2)=1;
        beta = reshape(V,sz(t2));
        
    else % if whole brain
        coords = [1:length(x0)]';
        V = df(fold).beta(:,coords,:);
        sz = size(V); t2 = sz~=1; t2(2)=1;
        beta = reshape(V,sz(t2));
%         beta  = squeeze(df(fold).beta(:,coords,:));
    end
    
    % recompute pred given params
    if size(pred_X_train,3) == 2
        pred_s = pred_X_train(coords,:,1); pred_s = num2cell(pred_s,2);
        pred_t = pred_X_train(coords,:,2); pred_t = num2cell(pred_t,2);
        pred_train = [pred_s pred_t];
        
        if makeTestPred == 1
            pred_s = pred_X_test(:,:,1); pred_s = num2cell(pred_s,2);
            pred_t = pred_X_test(:,:,2); pred_t = num2cell(pred_t,2);
            pred_test = [pred_s pred_t];
        end
        
    else
        pred_s = pred_X_train(coords,:,1); pred_s = num2cell(pred_s,2);
        pred_train = pred_s;
        
        if makeTestPred == 1
            pred_s = pred_X_test(:,:,1); pred_s = num2cell(pred_s,2);
            pred_test = pred_s;
        end
    end
    
    predictors_train = []; predictors_test=[];
    for vx = 1:size(pred_train,1)
        tmp1 = (cell2mat(pred_train(vx,:)')');
        predictors_train{vx,1} = [tmp1 df(fold).train_trend(:,df(fold).train_dcid)];
%         predictors_train{vx,1} = [cell2mat(pred_train(vx,:)')'];

        if makeTestPred == 1
            tmp2= cell2mat(pred_test(vx,:)')';
            predictors_test{vx,1}  = [tmp2  df(fold).test_trend(:,df(fold).test_dcid)];

%             predictors_test{vx,1}  = (cell2mat(pred_test(vx,:)')');

        end
        
    end
    
%     df(fold).train_trend(:,df(fold).train_dcid))
%     bb = num2cell(beta,2);
    
    if size(pred_X_train,3) == 2
        bb1 = beta(:,[1:2,df(fold).train_dcid+2]);
        bb2 = beta(:,[1:2,df(fold).test_dcid+2]);

%         bb = beta(:,1:2);

    else
         bb1 = beta(:,[1,df(fold).train_dcid+1]);
         bb2 = beta(:,[1,df(fold).test_dcid+1]);

%         bb = beta(:,1);
    end
    bb1 = num2cell(bb1,2);
    bb1 = cellfun(@transpose,bb1,'UniformOutput',false);

    trainPred = cell2mat(cellfun(@(x,y) x*y, predictors_train, bb1, 'UniformOutput',false)');
    varexp = calculateR2(trainData_raw',trainPred');

    % up date df
    df(fold).train_pred = trainPred;
    df(fold).varexp = varexp;
    
    if makeTestPred == 1
        bb2 = num2cell(bb2,2);
        bb2 = cellfun(@transpose,bb2,'UniformOutput',false);
        testPred  = cell2mat(cellfun(@(x,y) x*y, predictors_test , bb2, 'UniformOutput',false)');
        cv_varexp = calculateR2(testData_raw',testPred');
        df(fold).test_pred  = testPred;
        df(fold).cv_varexp = cv_varexp;
    end
    
end

end