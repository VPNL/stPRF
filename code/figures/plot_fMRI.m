clear; clc;


stPRFdir = pwd; % set path...
addpath(genpath('./code'))
download_toolboxes('toolbox_urls.txt', 'code/toolbox')

targetFile =   getAllFiles2('./results/paperfMRI/tables','*.mat',1);
plotDir = fullfile('./results/paperfMRI/figures');
if ~isfolder(plotDir); mkdir(plotDir); end

% select ROIs ====================================
roiList = {'V1','V2','V3','hV4','VO','LO','TO','V3AB','IPS'};

% extra ====================================
targetTmodel =  {'1ch_glm','1ch_dcts', '3ch_stLN'};

saveFig = 0;

%% VarExp
load(targetFile{3}) % load cross validated data

params=[];
[fn,saveName,stats] =  st_plot_varexp_violin(DT, params,roiList);
if saveFig == 1
    printnice(fn,[1 300],plotDir,saveName);
    printnice(fn,0,plotDir,saveName);
end
clear DT;

%% HRFs
load(targetFile{4}) % load hrftable
[fn,saveName] = st_plot_optHRF(hrfTable,sumTable);
if saveFig == 1
    printnice(fn,[1 300],plotDir,saveName.fig1);
end
clear hrfTable sumTable;
%% scatter - compare models
load(targetFile{2}) % load 9-run-concat optimized HRF data
df= defaultThresh(DT,0);
df = grabCommonVoxels(df,targetTmodel);
params=[]; 
for targetVar =  ["phase","sigma","ecc"]
    if targetVar =="sigma"
        [fn,saveName] =  st_plot_param_scatter_sigma(df, params,roiList, convertStringsToChars(targetVar));
    else
        [fn,saveName] =  st_plot_param_scatter(df, params,roiList, convertStringsToChars(targetVar));
    end
    if saveFig == 1
        printnice(fn,[1 300],plotDir,saveName)
        printnice(fn,0,plotDir,saveName)
    end
end
%% temporal param plot -- t2p
modelNumber=3;
params=[];
df= defaultThresh(DT,0);
[fn,saveName]=st_plot_temporal_bar(df, params,roiList,modelNumber);
if saveFig == 1
    printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
    printnice(fn.fig2,[1 300],plotDir,saveName.fig2);
    printnice(fn.fig3,[1 300],plotDir,saveName.fig3);
    printnice(fn.fig4,[1 300],plotDir,saveName.fig4);
    printnice(fn.fig5,[1 300],plotDir,saveName.fig5);
end

%% reconstruct and plot spatiotemporal profile --- all
% create spatiotemporal profile -> this is precomputed and uploaded on osf
% df= defaultThresh(DT,0); 
% channels=  {'sustained','transient'};
% for i = 1:length(channels)
%     channel = channels{i};
%     st_recon_stprf(df,params,roiList,channel,plotDir)
% end

df= defaultThresh(DT,0);
channels=  {'sustained','transient'};
for i = 1:length(channels)
    channel = channels{i};
    st_recon_plot(df,[],roiList,channel,plotDir)
end

%% reconstruct and plot spatiotemporal profile --- fov vs peri
% create spatiotemporal profile -> this is precomputed and uploaded on osf
% band = [0 4 ; 4  12];
% for b =1:2
%     df= defaultThresh(DT,0);
%     thresh = [];
%     thresh.ecc = band(b,1);  % min ecc
%     df = tableThreshold(df,thresh,'lower');
%     thresh = [];
%     thresh.ecc = band(b,2);  % max ecc
%     df = tableThreshold(df,thresh,'upper');
%     if b == 1
%         st_recon_stprf_ecc(df,[],roiList,'fov',plotDir);
%     elseif b ==2
%         st_recon_stprf_ecc(df,[],roiList,'peri',plotDir);
%     end
% end

% plot it (fov vs peri)
[fn,saveName] = st_recon_plot_ecc(DT,[],roiList,plotDir);
if saveFig == 1
    printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
    printnice(fn.fig1,0,plotDir,saveName.fig1);
end

%% temporal window vs spatial prf size
targetModel=3;
df= defaultThresh(DT,0);
[fn,saveName] = st_plot_sigma_window(df,[],roiList,targetModel,0);
if saveFig == 1
    printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
    printnice(fn.fig1,0,plotDir,saveName.fig1);
    printnice(fn.fig2,[1 300],plotDir,saveName.fig2);
    printnice(fn.fig2,0,plotDir,saveName.fig2);
end

%% ecc vs temporal window
targetModel = 3;
df= defaultThresh(DT,0);
testrois = {'V1','V2','V3'};
[fn,saveName,lme_temporal] = st_plot_temporal_ecc(df, [], testrois ,targetModel);
if saveFig == 1
    printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
end

%% sustained, transient contribution
df = defaultThresh(DT,1);
[fn,saveName]=st_plot_beta(df, params,roiList);
if saveFig == 1
    printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
    printnice(fn.fig2,[1 300],plotDir,saveName.fig2);
    printnice(fn.fig3,[1 300],plotDir,saveName.fig3);
end

%% ST-compare different measurements/methods
targetModel=3;
params=[];
df= defaultThresh(DT,0);
st_compareMethods(df, params,roiList,targetModel);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotExtra = 0;
if plotExtra
    %%  [extra] temporal param plot -- window
    df= defaultThresh(DT,0);
    [fn,saveName]=st_plot_temporal_bar(df, [],roiList,3,'window');
    
end