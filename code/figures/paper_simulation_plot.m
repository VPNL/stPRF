clear; clc;

simDir = pwd; % set path...
addpath(genpath('./code'))

% download required toolboxes
download_toolboxes('toolbox_urls.txt', 'code/toolbox')
%% download results file from OSF
% https://osf.io/3gwhz/files/osfstorage

%% set project name
simName = 'paperSim'; 
resultsDir = fullfile('results',simName);
sessionDir = fullfile(resultsDir,'data');

%% plot
% uses old model naming conventions
% 3ch_stLN == CST
% 1ch_glm == spatial
% 1ch_dcts == DN-ST

cd(fullfile(simDir,sessionDir))
estimated_table = 'dataTable/simulation/gray/df_paper.mat';
pt = load(estimated_table); % predicted table
gt = load(fullfile(simDir,sessionDir,'GT.mat')); % ground truth table
st_simulation_plot(gt,pt)
