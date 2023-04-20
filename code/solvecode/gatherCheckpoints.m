function gatherCheckpoints(vw,params,wProcess)
%%
checkdate = Constants.getDir.sessionDate;
tm = params.analysis.temporal.model;
tm = strrep(tm,'-','_');
checkDir  = getAllFiles(sprintf('./checkpoint/%s',checkdate),...
    'checkpoint.mat',1);
checkDir = checkDir(contains(checkDir,tm));
checkDir = checkDir(~contains(checkDir,'all'));

allCords = []; allParams=[];

if isempty(checkDir)
   return; 
end

savecheckDir = fullfile('./checkpoint',params.matFileName{end});
if ~isfolder(savecheckDir); mkdir(savecheckDir); end

if isfile(fullfile(savecheckDir,'checkpoint.mat'))
    return;
end



%%
for cf = 1:length(checkDir)
    coordsIndex =[];
    lcf = load(checkDir{cf});
    [~,roiName]=fileparts(fileparts(checkDir{cf}));
    roiName = strsplit(roiName,'-');
    roiName = roiName(contains(roiName,{'lh','rh','RH,LH'}));
    roiName = roiName{1};
    
    roiFiles = getAllFiles('./3DAnatomy/ROIs/' ,'*_toon.mat',2);
    if length(roiFiles(contains(roiFiles,roiName)))  == 1
        vw=loadROI(vw,roiFiles(contains(roiFiles,roiName)),[],[],1,0);
        coords = vw.ROIs(cf).coords;
        [coordsIndex, coords] = roiIndices(vw, coords);
        
    elseif length(roiFiles(contains(roiFiles,roiName))) ~= 1
        error('check_roi_name')
    end
    
    if size(lcf.estimatedParams,2) == numel(coordsIndex)
        coordsIndex = coordsIndex(~isnan(lcf.estimatedParams(1,:)));
        estimatedParams = lcf.estimatedParams(:,~isnan(lcf.estimatedParams(1,:)));
    else
        error('computed number of coords does not match')
    end
    
    allCords{cf}  = coordsIndex;
    allParams{cf} = estimatedParams;
end

%% gather and save
coordsIndex = cell2mat(allCords');
estimatedParams = cell2mat(allParams);

estimatedParams = NaN(size(estimatedParams,1),numel(wProcess));
estimatedParams(:,coordsIndex) = cell2mat(allParams);
save(fullfile(savecheckDir,'checkpoint.mat'),'estimatedParams','wProcess')


% wProcess = setdiff(wProcess,coordsIndex);


end

