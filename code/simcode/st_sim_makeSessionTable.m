function saveName = st_sim_makeSessionTable(sessionDir,rmFiles,viewType)
%% create table for fMRI pRF
if notDefined('sessionDir')
    sessionDir = pwd;
end

if notDefined('viewType')
    viewType = 'gray';
end

if strcmp(viewType, 'inplane')
    viewType = 'Inplane';
elseif  strcmp(viewType, 'gray')
    viewType = 'Gray';

end
% viewType = lower(viewType);
%
% if notDefined('makePred')
%     makePred = 0;
% end

%
% hg = initHiddenGray();
%
% roi = tc_roiStruct(hg, roi);
% [pred, RF, rfParams, variance_explained, blanks] = rmPlotGUI_makePrediction(M, coords);
% [M.tSeries, M.coords, M.params] = rmLoadTSeries(hg, params, roi, preserveCoords);

%%
cd(sessionDir)

dateStamp = Constants.getDir.sessionDate;
disp(['processing analysis file on:  ' dateStamp])

% rmFiles = getAllFiles(strcat('./',viewType,'/MotionComp_RefScan1/',dateStamp,'/'),'*sFit.mat',1);
[subpath,~] = fileparts(sessionDir);
[~,subjID] = fileparts(subpath);

mkdir('dataTable')
mkdir(['dataTable/' dateStamp '/' lower(viewType)])

% rmFiles = getAllFiles([stRootPath '/results/no'] ,'*st_*Fit.mat',2);
% rmFiles = getAllFiles([stRootPath '/results/no'] ,'*240*2ch*st*.mat',2);

% rmFiles = getAllFiles([stRootPath '/results/low'] ,'*2ch-exp-sig*st*.mat',2);
%%
for i = 1 :length(rmFiles)
    df2=[]; 
    
    rmFile = rmFiles{i};
    [~, b]=fileparts(rmFiles{i});
    fprintf('\n[st_makeSessionTable] [%d]: %s\n',i,b);
    tableName = fullfile('dataTable',dateStamp,lower(viewType),[b '-' num2str(3) '-df.mat']);
    if isfile(tableName); continue; end
    %     tableName = '';
    
    %     load(rmFile, 'model');
    load(rmFile, 'df');
    %%
    % get defualt params
    numFolds = df(1).numFolds;
    for fold = 1:numFolds
        
        tableName = fullfile('dataTable',dateStamp,lower(viewType),[b '-' num2str(fold) '-df.mat']);
        if isfile(tableName)
            fprintf('\n[st_makeSessionTable] [%d]: %d skipping\n',i,fold);
            continue;
        else
            fprintf('\n[st_makeSessionTable] [%d]: %d\n',i,fold);
        end
        % Naming
        modelname = convertCharsToStrings(b);
        pname = strsplit(modelname,'-');
        p.model  = pname{1};
        p.tmodel = pname{2};
        p.cond   = pname{3};
        p.roiName = pname{5};
        p.stage   = pname{end};
        p.split = fold;
        
        if contains(p.tmodel,'2ch') || contains(p.tmodel,'3ch')
            nchan = 2;
        else
            nchan = 1;
        end
        
            
        % if ROI is given or it is from an inplane view
        if ~strcmp(p.roiName,'all') && ~strcmpi(viewType,'Inplane')
            coords = df(fold).coordsIndex;
        else % if whole brain
            coords = [1:length(df(fold).x0)]';
        end
        
        if isempty(coords)  && ~strcmp(p.roiName,'all')
            coords = find(df(fold).searchFit)';
            if isempty(coords)
                error('check your coords')
            end

%             roiFiles = getAllFiles('./3DAnatomy/ROIs' ,sprintf('*%s.mat',p.roiName ),1);
%             if length(roiFiles) ~=1
%                 error('check your ROIs')
%             end
%             vw = initHiddenGray;
%             vw=loadROI(vw,roiFiles,[],[],1,0);
%             [coords, ~] = roiIndices(vw, vw.ROIs(1).coords);                
        end
            
        df2{i}.coords = single(coords);

        % deal with global fields
        fields = {'x0','y0','sigma','searchFit','varexp','cv_varexp','exponent'};
%         fields = {'searchFit'};
%         fields = {'varexp','cv_varexp'};

        for ff =1:length(fields)
            if isfield(df,fields{ff})
                each_param = df(fold).(fields{ff})';
                
                % only get major sigma cause we do not model DoG
                if strcmp(fields{ff},'sigma')
                    each_param = each_param.major';
                end
                
                % it is hack for uneven numbers of searchFit
                % need to fix it later...
                if  length(each_param) ~= length(df(fold).x0)
                    each_param = [each_param; zeros(length(df(fold).x0) - length(each_param),1)];
                end
                
            elseif strcmp(fields{ff},'searchFit')
                each_param = 0;
                each_param = repelem(each_param,length(df(fold).x0))';
            else
                error("unkownparam")
            end
            
          
            % restrict it to according coordinates
            each_param   =  each_param(coords);
            df2{i}.(fields{ff})   = single(each_param);
        end
        
        % get temporal fields
        switch lower(p.tmodel)
            case {'glm','1ch-glm','1ch_glm','spatial'}
                temproral_fields ={};
            case {'dn','1ch-dcts','1ch_dcts','dn_st'}
                temproral_fields = {'tau1', 'weight', 'tau2', 'nn', 'delay'};
            case {'2ch','2ch-exp-sig','2ch_exp_sig'}
                temproral_fields = {'tau_s', 'tau_ae', 'Lp', 'Kp', 'Kn'};
            case {'3ch','3ch_stln','cst'}
                temproral_fields = {'tau_s', 'tau_t'};
        end
        
        temporal=[];
        for ff =1:length(temproral_fields)
            if isfield(df,temproral_fields{ff})
                each_param = df(fold).(temproral_fields{ff})';
            else
                % if it is not in the dataField get default temporal params
                [prm,temporal_field]=getDefualtTemporalParams(p.tmodel);
                each_param = prm(strcmp(temporal_field,temproral_fields{ff}));
                each_param = repelem(each_param,length(df(fold).x0))';           
            end
            
            % restrict it to according coordinates
            each_param   =  each_param(coords);
            temporal.(temproral_fields{ff})   = single(each_param);
        end
        
        if isempty(temproral_fields)
            each_param = zeros(length(coords),1);
            temporal.shift =  each_param; % put in filler for the linear model;
        end

        df2{i}.temporal =  num2cell(cell2mat(struct2cell(temporal)'),2);
        
      

        
        % get beta values 
        beta  = single(squeeze(df(fold).beta(:,coords,1:nchan)));
        if nchan == 1
            df2{i}.beta =num2cell(beta',2);
        elseif nchan==2
            df2{i}.beta =num2cell(beta,2);
        end
        
        % get tc data
        trainData_raw = df(fold).train_data_raw';
        trainData_raw = trainData_raw(coords,:);
        testData_raw = df(fold).test_data_raw';
        testData_raw = testData_raw(coords,:);
        
        trainpred = df(fold).train_pred';
        trainpred = trainpred(coords,:);

        testpred = df(fold).test_pred';
        testpred = testpred(coords,:);


        df2{i}.trainData = num2cell(single(trainData_raw),2);
        df2{i}.testData = num2cell(single(testData_raw),2);
        df2{i}.trainPred = num2cell(single(trainpred),2);
        df2{i}.testPred  = num2cell(single(testpred) ,2);
        
      
        % turn it into a table
        T = struct2table(df2{i});
        
        % add string related variables
        if isrow(df(fold).train_set)
            df(fold).train_set = df(fold).train_set';
            df(fold).test_set = df(fold).test_set';
        end
        T.trainSet(:,1) = num2cell(df(fold).train_set,1);
        T.testSet(:,1) = num2cell(df(fold).test_set,1);
        
        T.model(:,1)  = convertCharsToStrings(p.model);
        T.tmodel(:,1) = convertCharsToStrings(p.tmodel);
        T.cond(:,1) = convertCharsToStrings(p.cond);
        T.roiName(:,1) = convertCharsToStrings(p.roiName);
        T.subjID(:,1) = convertCharsToStrings(subjID);
        T.split(:,1) = convertCharsToStrings(p.split);
        T.stage(:,1) = convertCharsToStrings(p.stage);

        save(tableName,'T','-v7.3');
        
        
        %     DT{i} = T;
        
        clear temp testData testData_a testData_b testData_c testPred
        clear testPred_a  testPred_b testPred_c pred pred_s pred_X
        clear predictors_test  predictors_test_a predictors_test_b predictors_test_c
    end
    
end

%% [search fit] model params only  
% datafiles = getAllFiles(fullfile('./dataTable/',dateStamp,'/',lower(viewType)) ,'*sFit*.mat',1);
datafiles = getAllFiles2(fullfile('./dataTable/',dateStamp,'/',lower(viewType)) ,'*sFit*.mat',1);

% fileList = getAllFiles2(dirName, fileExtension, appendFullPath)
for i = 1:size(datafiles,1)
    disp(['Table concating:   ' dateStamp '   -----  ' 'fileNumber:  '   num2str(i)])
    

    T1 = load(datafiles{i});
    
    % remove TC data for this
%     T1.T= removevars(T1.T,{'trainData','testData','trainPred','testPred','trainSet','testSet'});
%     T1.T= removevars(T1.T,{'trainData','testData','trainPred','testPred','trainSet','testSet'});

    DTavg{i} = T1.T;    

end
DT = vertcat(DTavg{:});
saveDir =['dataTable/' dateStamp '/' lower(viewType) '/'];
saveName= [saveDir 'df_' subjID];
save(saveName,'DT','-v7.3');

% %% [search fit] entire dataset including Timecourse dataset
% datafiles = getAllFiles(fullfile('./dataTable/',dateStamp,'/',lower(viewType)) ,'*sFit*.mat',1);
% DTavg = [];
% for i = 1:size(datafiles,1)
%     disp(['Table concating:   ' dateStamp '   -----  ' 'fileNumber:  '   num2str(i)])
%     
% 
%     T1 = load(datafiles{i});
%     DTavg{i} = T1.T;    
% 
% end
% DT = vertcat(DTavg{:});
% saveDir = ['../../../results/tables/',dateStamp,'/' ,lower(viewType) ,'/'];
% mkdir(saveDir)
% save([saveDir 'df_all' subjID],'DT','-v7.3');
    %%
%%%%%%

% below stuff are when you need to clean cell -> array
% or when you need to avg across splits.
% I decided to do this later when plotting. Just to have everything in one file

%%%%%%
% %
% % %% avg - compile
% %
% % % datafiles = getAllFiles('./dataTable/' ,'*avg*.mat',1);
% % datafiles = getAllFiles('./dataTable/' ,'*all*.mat',1);
% %
% %
% % TC=[]; DTavg=[];
% % for i = 1:size(datafiles,1)
% %     i
% %
% %     T1 = load(datafiles{i});
% %     T = T1.T;
% %     newT=[];
% %
% %     % unpack cells to arrays to mean them
% %     targetVariables = {'coords','x0','y0','sigma','varexp','cv_varexp','beta', ...
% %          'trainPred','testPred','trainData','testData'};
% % %     targetVariables = {'coords','x0','y0','sigma','varexp','beta', ...
% % %     'trainPred';
% %
% %     for vn = T.Properties.VariableNames
% %         paramName=vn;
% %         if max(ismember(targetVariables,paramName))
% %             if iscell(T.(paramName{:}))
% %                 newT.(paramName{:}) = cell2mat(T.(paramName{:}));
% %             else
% %                 newT.(paramName{:}) = T.(paramName{:});
% %             end
% %         end
% %
% %     end
% %     newT=struct2table(newT);
% %
% %     % Fix Beta for
% %     newT.beta = num2cell(newT.beta,2);
% %
% %     % add some important variables back to the table
% %     targetVariables = {'model','tmodel','cond','shuffle','roiName','subjID'};
% %     for vn = T.Properties.VariableNames
% %         paramName=vn;
% %         if max(ismember(targetVariables,paramName))
% %             newT.(paramName{:}) = T.(paramName{:})(1:height(newT));
% %         end
% %
% %     end
% %     %
% %     %
% %     DTavg{i} = newT;
% %
% %     tc=[];
% %     tc.name = 'fold3';
% %     targetVariables = {'trainSet','testSet'};
% %     for vn = 1:length(targetVariables)
% %         paramName=targetVariables(vn);
% %         t = [T.(paramName{:}){1}];
% %         tc.(paramName{:}) = t;
% %     end
% %
% %     TC{1} = tc;
% %
% %
% %
% % end
% %
% %
% % % DT = vertcat(DTavg{:});
% % DT = vertcat(DTavg{:});
% %
% % saveDir = '../../../results/tables/';
% % save([saveDir 'df_avg_' subjID],'DT','TC','-v7.3');
% %
% % %% compress splits and average run-folds
% %
% % datafiles = getAllFiles('./dataTable/' ,'*all*.mat',1);
% % % avgIdx = strfind(datafiles,'avg');
% % % avgIdx = find(~cellfun(@isempty,avgIdx));
% % % datafiles(avgIdx) = [];
% %
% % uniqueModels = cellfun(@(s) strsplit(s,'-'), ...
% %     getAllFiles2('./dataTable/' ,'*all*.mat',1), 'UniformOutput', false);
% % uniqueModels = cellfun(@(x) x(:,2), ...
% %     uniqueModels, 'UniformOutput', false);
% % uniqueModels = unique(cell2table(uniqueModels));
% %
% % nSplit = length(datafiles)/size(uniqueModels,1);
% % split_group = 1:length(datafiles);
% % split_group = reshape(split_group,nSplit,[])';
% %
% % for i = 1:size(split_group,1)
% %     i
% %     [s1,s2,s3]=datafiles{split_group(i,:)};
% %
% %     T1 = load(s1); T1 = T1.T;
% %     T2 = load(s2); T2 = T2.T;
% %     T3 = load(s3); T3 = T3.T;
% %
% %     T = [T1;T2;T3];
% %
% %     % unpack cells to arrays to mean them
% %     targetVariables = {'varexp','cv_varexp','cv_varexp_a','cv_varexp_b','cv_varexp_c','beta',...
% %         'trainPred','testPred','trainData','testData'};
% %     for vn = T.Properties.VariableNames
% %         paramName=vn;
% %         if max(ismember(targetVariables,paramName))
% %             if iscell(T.(paramName{:}))
% %                 T.(paramName{:}) = cell2mat(T.(paramName{:}));
% %             else
% %                 T.(paramName{:}) = T.(paramName{:});
% %             end
% %         end
% %
% %     end
% %
% %     DT{i} = T;
% %
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%% mean %%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % mean across splits!!
% %     %         meanT=[];
% %     %         meanT= varfun(@mean,T,'InputVariables',...
% %     %             {'x0','y0','sigma','varexp','beta'} ,...
% %     %             'GroupingVariables','coords');
% %     % %                 meanT= varfun(@mean,T,'InputVariables',...
% %     % %             {'x0','y0','sigma','varexp','cv_varexp','beta'} ,...
% %     % %             'GroupingVariables','coords');
% %     %
% %     %         % clean up, remove 'mean_' from names
% %     %         for cn = 1:length(meanT.Properties.VariableNames)
% %     %             nn = meanT.Properties.VariableNames{cn};
% %     %             meanT.Properties.VariableNames{cn} = erase( nn , 'mean_' );
% %     %         end
% %     %
% %     %         % put beta back into cell format and remove groupcount header
% %     %         meanT.beta = num2cell(meanT.beta,2);
% %     %
% %     %         meanT.GroupCount = [];
% %     %
% %     %         % add some important variables back to the table
% %     %         targetVariables = {'coords','model','tmodel','cond','shuffle','roiName','subjID' };
% %     %
% %     %         for vn = T.Properties.VariableNames
% %     %             paramName=vn;
% %     %             if max(ismember(targetVariables,paramName))
% %     %                 meanT.(paramName{:}) = T.(paramName{:})(1:height(meanT));
% %     %             end
% %     %
% %     %         end
% %     %
% %     % %            'trainPred','testPred','trainData','testData'
% %     %
% %
% %
% %
% %
% %
% % end
% %
% %
% % DT = vertcat(DT{:});
% % saveDir = '../../../results/tables/';
% % save([saveDir 'df_' subjID],'DT','-v7.3');

%% tc-related table
%
% datafiles = getAllFiles('./dataTable/' ,'*all*.mat',1);
%
%
% if length(datafiles)==18
%     split_group = 1:18; % hard coded
%     split_group = reshape(split_group,3,6)';
% end
%
% if length(datafiles)==18
%     for i = 1:size(split_group,1)
%         [s1,s2,s3]=datafiles{split_group(i,:)};
%
%         T1 = load(s1); T1 = T1.T;
%         T2 = load(s2); T2 = T2.T;
%         T3 = load(s3); T3 = T3.T;
%
%         T = [T1;T2;T3];
%
%         i
%         tc = [];
%         targetVariables = {'model','tmodel','cond','shuffle','roiName','subjID'};
%         for vn = T.Properties.VariableNames
%             paramName=vn;
%             if max(ismember(targetVariables,paramName))
%                 tc.(paramName{:}) = T.(paramName{:})(1:height(T1));
%             end
%         end
%         tc=struct2table(tc);
%
%
%         targetVariables = {'trainData','testData','trainPred','testPred'};
%         for vn = 1:length(targetVariables)
%             paramName=targetVariables(vn);
%             t = [T1.(paramName{:}) T2.(paramName{:}) T3.(paramName{:})];
%             t = cellfun(@transpose,t,'UniformOutput',false);
%             t = cellfun(@cell2mat,num2cell(t,2),'un',0);
%             tc.(paramName{:}) = t;
%         end
%           TC{i} = tc;
%     end
% end
%
% %
% else
%     for i =1 :length(datafiles)
%         i
%         dfile = datafiles{i};
%         load(dfile);
%
%         split_group(1,:)
%
%         DT{i} = T;
%
%     end
%
% end

%%
%

% for cc = 1:size(DT,2)
%     a=table2array(DT(:,cc));
%     if iscell(a)
%         b=cellfun(@single,a,'un',0);
%         DT(:,cc)=array2table(b);
%     else
%       DT(:,cc)=  array2table(single(a));
%     end
% end

%%

%
% % a=DT(:, contains(DT.Properties.VariableNames, 'noise'))
%
% DT = vertcat(DT{:});
% save('T_low_240.mat','DT');
% % tblB = sortrows(T,{'varexp'},{'descend'})
% % a = DT(DT.dur == "240" & DT.RF == "1",:);
% a = DT(DT.solvestage == "gFit" | DT.solvestage == "sFit",:);


end