function st_prfRun(sessionDir,roifile,stimfile,analysisoption,prfModel,temporalModel,prefix,wSearch,userData,useGPU)
%
% Step 4 of the workflow for analyzing stPRF
% Runs the spatiotemporal prf model 
%
% Default Input valuesw

%
% KIS 2020 (adapted from code by RL & JG & DF)

%% TODO:
% remove: stimseq
% now everything is based on input stimulus design

%%
% is this Sherlock?
% if isfile('/scratch/users/insubkim/yesSherlock.txt')
%     tbUseProject('stRet_toolbox');
% else
% %     tbUseProject('stprf_toolbox');
% %     addpath(genpath('/Volumes/Samsung_T5/project/'));
% end


%% Default inputs


if notDefined('sessionDir')
    sessionDir= pwd;
end

if notDefined('roifile')
    doROI = false;
else
    doROI = true;
end

if notDefined('stimfile')
%     stimfile = fullfile('.', 'Stimuli', 'images_and_params');
    stimfile = getAllFiles('./Stimuli','images_and_params*',1);

end

if notDefined('prfModel') % {'one oval gaussian' | 'onegaussian' | 'css'}
    prfModel = {'st'};
else
    prfModel = cellstr(prfModel); 
end

if notDefined('temporalModel') % {'1ch-glm' | '1ch-dcts' |'2ch-exp-sig'};
    if strcmp(prfModel{1},'st')
        temporalModel = '2ch-exp-sig';
    else
        temporalModel = '1ch-glm';
    end
end

% 1 = crossValidation 
% 99 = Development
if notDefined('analysisoption')
    analysisoption = 1;     
end

if notDefined('prefix')
    prefix = 'subj01';
end

if notDefined('wSearch')
    % search type.
    % 1 = grid search only ("coarse"),
    % 2 = minimization search only ("fine"),
    % 3 = grid followed by minimization search [default]
    % wSearch = 3;    
    wSearch = 8;
end


if notDefined('userData')
    userData = [];
end


[~,cmdout] = system('echo $HOSTNAME');
%
if notDefined('useGPU')
    % isfile('/scratch/users/insubkim/yesSherlock.txt') ||
    if contains(cmdout,'new-nefesh') || contains(cmdout,'strelka')  || contains(cmdout,'henrygibbons')
        %    useparallel = 1;
        useGPU =  1;
    elseif ~isempty(regexp(cmdout,'sh([0-9])+')) && regexp(cmdout,'sh([0-9])+') == 1
        fprintf('+*+*+*+*+**+*+*+*+*+*+*+*+ GPU ON \n +*+*+*+*+**+*+*+*+*+*+*+*+')
        useGPU =  1;
        fprintf('+*+*+*+*+**+*+*+*+*+*+*+*+ GPU ON \n +*+*+*+*+**+*+*+*+*+*+*+*+')
    else
        useGPU =  0;
    end
end

% if notDefined('hrfver')
%     hrfver =[];
% end


useparallel= 0;
if useparallel
    % open parpool with 12 workers
    openParPool(12)
end

%% Setup  ----------------------
close all; clc;

cd(sessionDir);

% load('ExpDesign.mat')
% T = cell2table( ExpDesign ) ;
% [C, ~, ~] = unique(T);
% for ec = 1:height(C)
%     list_rmName{ec} = strcat(C.ExpDesign1{ec}, '_',C.ExpDesign2{ec});
% end

% [~, sessioncode] = fileparts(sessionDir);

%% key attribute settings  ----------------------
% radius of circle retinotopy in visual angle degrees
stimradius = 12;

% detrending frequeny maximum (cycles per scan):
% 1 means 3 detrending terms, DC (0 cps), 0.5, and 1 cps
detrend = 1;

%%
optionBundle = st_selectOption(analysisoption);

if strcmp(prfModel,'st')
    stimType = 'StiminMS';
else
    stimType = 'StimFromScan';
end


%% Run the prf model  ----------------------

if isfile('/scratch/users/insubkim/yesSherlock.txt')
    anatPath = [sessionDir '/3DAnatomy' ];
    setVAnatomyPath([anatPath,'/t1.nii.gz'])
end
% open the session
vw = initHiddenGray;
% vw = initHiddenInplane;

% need some global variables later
load mrSESSION;

%% loop over the datatypes  ----------------------

% set current dataTYPE
vw = viewSet(vw, 'curdt', optionBundle.rmName);

% get the dataType struct
dtstruct = viewGet(vw, 'dtstruct');

% get the data type number for later
dataNum = viewGet(vw, 'curdt');

if doROI == 1
    vw=loadROI(vw,fullfile('./3DAnatomy','ROIs',roifile),[],[],1,0);
%   loadROI(vw, filename, select, clr, absPathFlag, local)
    fprintf('\n[st_prfRun] loading ROI: %s\n',roifile)
end


%% create stimulus  ----------------------

%%%% if NIFTI is given as stimulus
%     [~, f, e] = fileparts(stimfile);
%     ni = niftiRead(fullfile(homedir, 'Stimuli', sprintf('%s%s', f, e)));
%
%     images = squeeze(ni.data);
%     images=single(images);
%
%     pixdim = niftiGet(ni, 'pixdim');
%     tr     = pixdim(end);
%     sprintf('/n/n USING TR:%2.2f/n/n',tr)
%%%%
% 
% mkdir(fullfile(sessionDir, 'Stimuli'));
% copyfile(stimfile, fullfile(sessionDir, 'Stimuli'));
% 
% load(stimfile);
% images = stim;
% stimulus.seq = 1:size(images,3);
% stimulus.seqtiming = (stimulus.seq-1);
% 
% if contains(rmNameSplit{2},'seq')
%     stimfileMat = fullfile('.', 'Stimuli', 'images_and_params');
%     save(stimfileMat, 'images', 'stimulus');
% else
%     stimfileMat = fullfile('.', 'Stimuli', 'images_and_params_shuffle');
%     save(stimfileMat, 'images', 'stimulus');
% end


%% set PRF params  ----------------------

params = rmCreateStim(vw);

for pp = 1 :length(params)
    params(pp).nFrames = viewGet(vw, 'nFrames');
    params(pp).framePeriod = viewGet(vw, 'framePeriod');
    tem.totalFrames = mrSESSION.functionals(1).totalFrames;
    % prescanDuration is going to be cliped
    params(pp).prescanDuration = (tem.totalFrames - params(pp).nFrames) * params(pp).framePeriod;
    params(pp).stimType = stimType;
    params(pp).stimSize = stimradius; % stimulus radius (deg visual angle)
    params(pp).nDCT     = detrend;
%     params(pp).split = wSplit;
    params(pp).imFile     = stimfile{pp};    % file containing stimulus images
    params(pp).paramsFile = stimfile{pp};    % file containing stimulus parameters
    params(pp).hrfType = 'two gammas (SPM style)';
    params(pp).hrfparams = {[1.6800, 3, 2.0500], [5.4000, 5.2000, 10.8000, 7.3500, 0.3500]};
    params(pp).imFilter   = 'none';
    params(pp).shuffled = false;      % need to fix this later ...............
    params(pp).fliprotate = [0, 0, 0]; %[L/R flip up/down flip rotate degrees]
    params(pp).stimWidth = 90;
    params(pp).stimStart = 0;
    params(pp).stimDir = 0;
    params(pp).nCycles = 1;
    params(pp).nStimOnOff = 0;
    params(pp).nUniqueRep = 1;
end

% params.nFrames = viewGet(vw, 'nFrames');
% params.framePeriod = viewGet(vw, 'framePeriod');
% tem.totalFrames = mrSESSION.functionals(1).totalFrames;
% 
% % prescanDuration is going to be cliped
% params.prescanDuration = (tem.totalFrames - params.nFrames) * params.framePeriod;
% 
% params.stimType = 'StiminMS';
% params.stimSize = stimradius; % stimulus radius (deg visual angle)
% params.nDCT     = detrend;
% params.split = wSplit;
% params.imFile     = stimfileMat;    % file containing stimulus images
% params.paramsFile = stimfileMat;    % file containing stimulus parameters
% params.hrfType = 'two gammas (SPM style)';
% params.hrfParams = {[1.6800, 3, 2.0500], [5.4000, 5.2000, 10.8000, 7.3500, 0.3500]};
% params.imFilter   = 'none';
% params.shuffled = false;      % need to fix this later ...............
% params.fliprotate = [0, 0, 0]; %[L/R flip up/down flip rotate degrees]
% params.stimWidth = 90;
% params.stimStart = 0;
% params.stimDir = 0;
% params.nCycles = 1;
% params.nStimOnOff = 0;
% params.nUniqueRep = 1;
%% getting parameter values for prf model fit

% store it
dataTYPES(dataNum).retinotopyModelParams = params;
% dataTYPES = dtSet(dataTYPES, 'rm stim params', sParams);

% save it
saveSession;
load('mrSESSION.mat');

% store params in view struct
vw = viewSet(vw, 'rmParams', params);

if wSearch == 7
    saveTmodelName = strrep(temporalModel,'-','_');
    rmFilename = getAllFiles(['./Gray/' optionBundle.rmName] ,sprintf('%s*%s*gFit.mat',prfModel{1},saveTmodelName),1);
    vw = viewSet(vw, 'rmfile', rmFilename{1});
end
%% Put the rm params into the view structure

%     vw = rmLoadParameters(vw);
    % the function rmLoadParameters used to call both rmDefineParameters
    % and rmMakeStimulus. If we do it here so that we can give it arguments
    % outside of the default (eg previously, sigma major and minor would be
    % identical despite having prfModel = {'one oval gaussian'} when
    % specifying it as an argument in vw = rmMain(vw, ...)

    % store params in view struct
    vw = viewSet(vw, 'rmParams', params);

fprintf('\n[st_prfRun] This is homedir: %s\n',sessionDir)
fprintf('\n[st_prfRun] This is stimfile of the first: %s\n',stimfile{1})
fprintf('\n[st_prfRun] This is stimradius: %i\n',stimradius)
fprintf('\n[st_prfRun] This is temporalModel: %s\n',temporalModel)


%% RUN THE PRF!

% name the ret model - whole brain
try
    outputDir = fullfile('./Gray/',optionBundle.rmName,Constants.getDir.sessionDate,'/');
    
    mkdir(outputDir);
    
    rmName = strrep(optionBundle.rmName,'_','-');
    saveTmodelName = strrep(temporalModel,'-','_');
    if doROI
        [~,roiname] = fileparts(roifile);
        outFileName = [Constants.getDir.sessionDate '/' prfModel{1},'-',saveTmodelName,'-',rmName, '-' roiname];
    else
        outFileName = [Constants.getDir.sessionDate '/' prfModel{1},'-',saveTmodelName,'-',rmName, '-all'];
    end

catch
    outputDir = fullfile('./Gray/',optionBundle.rmName,'/');
    
    mkdir(outputDir);
    
    rmName = strrep(optionBundle.rmName,'_','-');
    saveTmodelName = strrep(temporalModel,'-','_');
    if doROI
        [~,roiname] = fileparts(roifile);
        outFileName = [outputDir '/' prfModel{1},'-',saveTmodelName,'-',rmName, '-' roiname];
    else
        outFileName = [outputDir '/' prfModel{1},'-',saveTmodelName,'-',rmName, '-all'];
    end

end

vw = rmMain(vw, [], wSearch, ...
    'model', prfModel, ...
    'matFileName', outFileName, ...
    'keepAllPoints', optionBundle.keepAllPoints, ...
    'numberStimulusGridPoints', optionBundle.numberStimulusGridPoints, ...
    'coarsetofine', optionBundle.coarsetofine, ...
    'coarsetblurparams', optionBundle.coarsetblurparams, ...
    'calcPC', optionBundle.calcPC, ...
    'prefix',prefix, ... 
    'cache', optionBundle.cache, ...
    'cv', optionBundle.cv, ...
    'useparallel', useparallel, ...
    'temporalModel',temporalModel, ...,
    'useGPU',useGPU, ...
    'userInputData', userData, ...
    'hrfver', optionBundle.hrfver);
%     'hrfver', optionBundle.hrfver);
%    'normrf', 1, ...
%     'normStimRF', optionBundle.normrf, ...

end
