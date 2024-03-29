%% spatiotemporal pRF Demo
% 1) download toolboxes
% 2) set variables and create JSON files
% 3) create Synthetic timecourses
% 4) Solve spatiotemporal PRF models
% 5) Plot results

%%
% download required toolboxes
download_toolboxes('toolbox_urls.txt', 'code/toolbox')

%% generate synthBOLD
clear; clc;
simDir = pwd; % set path...
addpath(genpath('./code'))

%% set variables and params

% set project name
simName = 'example'; noiselevel = 'noise1'; 
% simName = 'match'; noiselevel = 1; 

resultsDir = fullfile('results',simName);
sessionDir = fullfile(resultsDir,'data');

% set noise level and stimulus size (radius, in visual angle, deg)
nVoxels = 10; StimSize = 12;  temporalSampleRate = 100; 
nRuns = 2; % using 1 run to reduce compute time otherwise
% use 9 runs to include all the temporal conditions in the experiment 

% set global variables
outputDir=fullfile(simDir,sessionDir);
setConstants(outputDir,outputDir,'simulation');



%% craete jsonfile
cd(simDir)
mkdir(fullfile('results',simName,'synthBOLD'));

% create expNames and stimFiles to use
expNames = string([]);
stimFiles = string([]);
for runNumber =  1:nRuns
    expNames(runNumber) = "ST_run" + runNumber;
    stimFiles(runNumber) = "./Stimuli/images_and_params_run0" + runNumber + ".mat";
end

% create randomized spatial RF params [x,y,sigma]
RF = randomRFparam(StimSize,nVoxels);

% create randomized extra params for a specfic model [tau, n]
tParam1 = randomTemporalparam('CST',nVoxels); % CST model
tParam2 = randomTemporalparam('DN-ST',nVoxels); % DN model

% load pre-defined jsonfiles
jsonTemplate = fullfile('template','template.json');

% manipulate Json files with these randomly generated params
jsonDir = fullfile(resultsDir,'jsonfiles');  mkdir(jsonDir);
randjson(jsonDir,jsonTemplate,nVoxels,expNames,stimFiles, ...
    StimSize,temporalSampleRate,RF,tParam1,tParam2,noiselevel);

%% synthetic BOLD generation
jFiles2 = getAllFiles2(jsonDir, sprintf('*.json'),2);
for jj = 1:length(jFiles2)
    json = jFiles2{jj};
    [~,jsonname]=fileparts(json);
    outputDir = fullfile(resultsDir,'synthBOLD',jsonname);
    DTcalc = synthBOLDgenerator(json,outputDir);
end
%% gather and clean BOLD files
removeFirstBlankSignals = 10;  % 10 TRs
rfiles = getAllFiles(fullfile(resultsDir,'synthBOLD'),sprintf('*table.mat'),1);
st_cleanDT(rfiles,resultsDir,removeFirstBlankSignals);

%% visualize synthetic timecourse
t = load(fullfile(resultsDir,'data/GT.mat'));
voxelNumber = 1; 
st_viewTC(t.DT,voxelNumber)

%% solve pRF model
% [Warning]: A GPU is highly recommended for this task due to
% the significant compute time required. 
% Running the task on a CPU may result in extremely long processing times,
% We strongly advise using a GPU to ensure that the task completes in a timely manner.

% set mrVista files
dataTemplate = fullfile('./template','data_template');
copyfile(dataTemplate,fullfile(resultsDir,'data'))
copyfile('./Stimuli',fullfile(resultsDir,'data','Stimuli'))

load(fullfile(resultsDir,'data','mrSESSION_template.mat'))
mrSESSION.functionals.totalFrames = 220;  %%%  need to set this
dataTYPES.scanParams.nFrames = mrSESSION.functionals.totalFrames - 10; %%%  need to set this
sp = dataTYPES.scanParams;
f =  mrSESSION.functionals;

for i =  1:nRuns
    dataTYPES.scanParams(i) = sp;
    mrSESSION.functionals(i) = f;
end
save(fullfile(resultsDir,'data','mrSESSION.mat'),'dataTYPES','mrSESSION','vANATOMYPATH')

%%
cd(simDir)
roifile = [];
stimfile = getAllFiles('./Stimuli','images_and_params*',1);
analysisoption = 5; % all
prfModel = 'st';
subjNumber=1;
wSearch = 8;

temporalModels = {'CST','DN-ST','spatial'};

% actually start solving
for tm = 1:length(temporalModels)
    cd(simDir)
    temporalModel =  temporalModels{tm};
    userData = sprintf('bold_%s_voxel-all.mat',temporalModel);
    st_prfRun(sessionDir,roifile,stimfile,analysisoption, ...
    prfModel,temporalModel,sprintf('subj%02d', subjNumber),wSearch, ...
    userData)

end

%% create results table
cd(fullfile(simDir,sessionDir))
rfiles = getAllFiles('./Gray/simulation_result/simulation', ...
                     sprintf('*sFit.mat'),2);
estimated_table = st_sim_makeSessionTable(fullfile(simDir,sessionDir),rfiles,'gray');
%% plot
cd(fullfile(simDir,sessionDir))
estimated_table = Constants.getDir.sim_pt;
pt = load(estimated_table); % predicted table
gt = load(fullfile(simDir,sessionDir,'GT.mat')); % ground truth table
st_simulation_plot(gt,pt)

%% 1.97 0.47
