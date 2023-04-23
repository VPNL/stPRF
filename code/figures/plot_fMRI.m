clear; clc;

stPRFdir = pwd; % set path...
addpath(genpath('./code'))

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
for targetVar = ["phase","sigma","ecc"]
    [fn,saveName] =  st_plot_param_scatter(DT, params,roiList, convertStringsToChars(targetVar));
    if saveFig == 1
        printnice(fn,[1 300],plotDir,saveName)
        printnice(fn,0,plotDir,saveName)
    end
end

%% temporal param plot -- t2p
modelNumber=3;
params=[];
df= defaultThresh(DT,0);a
[fn,saveName]=st_plot_temporal_bar(df, params,roiList,modelNumber);
if saveFig == 1
    printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
    printnice(fn.fig2,[1 300],plotDir,saveName.fig2);
    printnice(fn.fig3,[1 300],plotDir,saveName.fig3);
    printnice(fn.fig4,[1 300],plotDir,saveName.fig4);
    printnice(fn.fig5,[1 300],plotDir,saveName.fig5);
end

%% reconstruct and plot spatiotemporal profile --- all
df= defaultThresh(DT,0);
channels=  {'sustained','transient'};
for i = 1:length(channels)
    channel = channels{i};
    st_recon_stprf(df,params,roiList,channel,plotDir)
    st_recon_plot(df,params,roiList,channel,plotDir)
end
%% reconstruct and plot spatiotemporal profile --- fov vs peri

% create spatiotemporal profile
% band = [0 4 ; 4  12];
% for b =2 
%     df= defaultThresh(DT,0);
%     thresh = [];
%     thresh.ecc = band(b,1);  % min ecc
%     df = tableThreshold(df,thresh,'lower');
%     thresh = [];
%     thresh.ecc = band(b,2);  % max ecc
%     df = tableThreshold(df,thresh,'upper');
%     if b == 1
%         st_recon_stprf_ecc(df,params,roiList,'fov',plotDir);
%     elseif b ==2
%         st_recon_stprf_ecc(df,params,roiList,'peri',plotDir);
%     end
% end

% plot it (fov vs peri)
st_recon_plot_ecc(DT,params,roiList,plotDir);


%% temporal window vs spatial prf size
params =[];
targetModel=3;
df= defaultThresh(DT,0);
[fn,saveName] = st_plot_sigma_window(df,params,roiList,targetModel);
if saveFig == 1
    printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
    printnice(fn.fig1,0,plotDir,saveName.fig1);
    printnice(fn.fig2,[1 300],plotDir,saveName.fig2);
    printnice(fn.fig2,0,plotDir,saveName.fig2);
end

%% ecc vs temporal window
params=[];
targetModel = 3;
df= defaultThresh(DT,0);
testrois = {'V1','V2','V3'};
[fn,saveName] = st_plot_temporal_ecc(df, params, testrois ,targetModel);
if saveFig == 1
    printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
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
    %% [extra] beta plot
    df = defaultThresh(DT,1);
    [fn,saveName]=st_plot_beta(df, params,roiList);
    if saveFig == 1
        printnice(fn.fig1,[1 300],plotDir,saveName.fig1);
        printnice(fn.fig2,[1 300],plotDir,saveName.fig2);
        printnice(fn.fig3,[1 300],plotDir,saveName.fig3);
    end
    
    %%  [extra] temporal param plot -- window
    modelNumber=3;
    params=[];
    df= defaultThresh(DT,0);
    [fn,saveName]=st_plot_temporal_bar(df, params,roiList,modelNumber,'window');
    
end