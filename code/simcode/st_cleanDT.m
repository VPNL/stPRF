function st_cleanDT(rfiles,resultsDir,removeInitial)
% clean synthBOLD into a nice looking table
% also reformat BOLD timecourse into to the mrVista format

% amount of initial signal (blank)
if exist('removeInitial','var')
else
    removeInitial = 0;
end

%% gather and Clean BOLD files

AllBOLD = []; AllBOLDpure=[]; GT =[];
% Regular expression pattern to match run numbers
pattern = 'run(\d+)_';
InputFileDelimit = '/';

% Extract run numbers from file paths using regular expression and convert to numeric array
runs = (cellfun(@(x) str2double(regexp(x, pattern, 'tokens', 'once')), rfiles))';
nRun = length(unique(runs));

% Loop over unique runs
for nr = 1:nRun
    
    voxelFile = selectFile(rfiles,sprintf('run%d',nr),InputFileDelimit);
    
    gt = []; DT =[]; voxel= [];
    for vf =  1:length(voxelFile)
        tmpDT=  load(voxelFile{vf});
        DT = [DT;tmpDT.DT] ;
        voxel = [voxel; repmat(vf,height(tmpDT.DT),1)];
    end
    
    
    % create groundtruth table
    SNR   = DT.SNR;
    Noise = DT.Noise.voxel;
    run   = repmat(nr,length(Noise),1);
    
    tcn = []; tc=[];
    for vx = 1: height(DT)
        tcn(vx,:) = DT(vx,:).pm.BOLDnoise;
        tc(vx,:) = DT(vx,:).pm.BOLD;
    end
    tcn = tcn(:,removeInitial+1:end);
    tc  = tc(:,removeInitial+1:end);
    
    T = table(voxel,SNR,Noise,run,tcn,tc);
    gt = [DT.RF(:,1:5) DT.Temporal T];
    GT = [GT; gt];
    
    [a, b]=fileparts(fileparts(voxelFile{1}));
    if exist(fullfile(fileparts(a),'data'), 'dir') ~= 7
        mkdir(fullfile(fileparts(a),'data'));
    end
    %     save(fullfile(fileparts(a),'voxels',strrep(b,'voxels001','voxels-all')),'DT','-v7.3');
    
    temporal = unique(DT.Temporal.temporalModel);
    for et = 1:length(temporal)
        Index = find(contains(DT.Temporal.temporalModel, temporal{et}));
        dt = DT(Index,:);
        for vx = 1:height(dt)
            AllBOLD(:,vx,nr,et) = dt(vx,:).pm.BOLDnoise(removeInitial+1:end);
            AllBOLDpure(:,vx,nr,et) = dt(vx,:).pm.BOLD(removeInitial+1:end);
        end
    end
end

for et = 1:length(temporal)
    BOLD = AllBOLD(:,:,:,et);
    
    % save mrVista BOLD
    saveName = sprintf('bold_%s_%s',temporal{et},['voxel-all']);
    save(fullfile(resultsDir,'data',saveName),'BOLD')
    sn = concatRuns(AllBOLD(:,:,:,et));
    s = concatRuns(AllBOLDpure(:,:,:,et));
    for i=1:size(sn,2)
        catSNR(i,et)  = snr(s(:,i), (sn(:,i) - s(:,i)));
    end
    
end



% clean & save Ground-truth datatable
% RF
[aGroup, voxel] = findgroups(GT.voxel);
RF = splitapply(@mean, [GT.Centerx0 GT.Centery0 GT.sigmaMajor], aGroup);
RF = table(voxel, RF);

% temporal params
[aGroup,  temporalModel] = findgroups(GT.temporalModel);
tp=[];
for eg = 1:length(temporalModel)
    out = cell(height(RF),1);
    [out{:}] = (GT(find(aGroup==eg),:).tParams.values);
    tp= [tp out];
end

% gather
DT = [];
for tm = 1:length(temporalModel)
    tparam = tp(:,tm);
    SNR    = catSNR(:,tm);
    tmodel = temporalModel(tm);
    tmodel = repmat(tmodel, length(tparam), 1);
    
    ntc = concatRuns(AllBOLD(:,:,:,tm))';
    tc = concatRuns(AllBOLDpure(:,:,:,tm))';
    
    T = table(tmodel,tparam,SNR,ntc,tc);
    DT{tm} = [RF T];
    
end

% save gound-truth Datatable
save(fullfile(resultsDir,'data',['GT']),'DT')

end