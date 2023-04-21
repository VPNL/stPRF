function [fn,saveName]=st_plot_sigma_window(DT,params, roiList,targetModel)



mc = getmycolors(3);

if targetModel == 2
    tms = {'\tau_1', '\itw', '\tau_2', '\itn_D_N', '\sigma_D_N'};
    nParam = 5;
elseif targetModel == 3
    tms = {'\tau', '\tau'};
    nParam = 1;
end

if notDefined('targetVar')
   targetVar = 'temporal';
end

params.modelNumber = targetModel;

% initial filter
df = tableFilter(DT,params);

% give number to ROIs
df = st_tablegiveROInumber(df, roiList);

%%
% threshold extreme end-points
thresh = [];
thresh.temporal = 99;
df = tableThreshold(df,thresh,'upper'); % remove anything upper;

thresh = [];
thresh.temporal = 4; 
df = tableThreshold(df,thresh,'lower'); % remove voxels above 12 ecc
%%
params=[];
params.modelNumber = targetModel;
df2 = tableFilter(df,params);

% % unpack table cells to arrays
df2 = tableArrayFormat(df2);
id = unique(DT.subjID);

for sj = 1:length(id)
    params.subjID =  id(sj);
    for eachROI = 1:length(roiList)
        params.roiNumber = eachROI;
        df3 = tableFilter(df2,params);

        x0{sj,eachROI} = df3.x0;
        y0{sj,eachROI} = df3.y0;
        sigma{sj,eachROI} = df3.sigma;
        ecc{sj,eachROI} = df3.ecc;
        window{sj,eachROI} = df3.window;
        expo{sj,eachROI} = df3.exponent;
        
        for tp = 1:nParam %size(df3.temporal,2)
            agg{sj,eachROI} = df3.(targetVar)(:,tp);
        end
    end
    
end

%%

%%
fn.fig1 = figure(1); clf;
saveName.fig1 = ['sigma_and_tau'];
% figure(1);
ScrSz = get(0, 'screensize');
set(gcf, 'position',  [ScrSz(1) ScrSz(2) ScrSz(3)/2 ScrSz(4)/1.5]);

for i =1:length(roiList)

    w = cellfun(@median, window(:,i));
    w =w(~isnan(w));
    x_e = nanstd(w) / sqrt(size(w, 1) - 1);
    x =  median(cell2mat(window(:,i)));
    
    s = cellfun(@median, sigma(:,i));
    s = s(~isnan(s));
    y_e = nanstd(s) / sqrt(size(s, 1) - 1);
    y =  median(cell2mat(sigma(:,i)));
    
    
    ev = cellfun(@median, expo(:,i));
    ev = ev(~isnan(ev));
    e_sem = nanstd(ev) / sqrt(size(ev, 1) - 1);
    ee =  median(cell2mat(expo(:,i)));
    
    errorbar(x,y,y_e,y_e,x_e,x_e,'o', ...
        'Color',mc(i, :),"MarkerFaceColor",mc(i, :),'markersize',(1/ee)*3,'linewidth',3); hold on;

    
end
xlim([4 15])
xticks([0:2:20]);
xticklabels(xticks*10);

tch_set_axes;
title('temporalWindow vs prfSize (SEM)')
% % legend(roiList,'orientation','horizontal','location','southoutside'); legend box off;
xlabel('temporal window (ms)');
ylabel('pRF size (deg)');
set(gca,'XMinorTick','on','yMinorTick','on');
% xlim([60 150])

%%

fn.fig2 =figure(2); clf;
saveName.fig2 = ['sigma_and_tau_CI'];
% figure(1);
ScrSz = get(0, 'screensize');
set(gcf, 'position',  [ScrSz(1) ScrSz(2) ScrSz(3)/2 ScrSz(4)/1.5]);

for i =1:length(roiList)
    
    
    s = (cell2mat(sigma(:,i)));
    w = (cell2mat(window(:,i)));
    
    bootstat = bootstrp(100,@median,[w s ]);
    x = mean(bootstat(:,1));
    y = mean(bootstat(:,2));
    
    alpha = 0.05;
    CI = prctile(bootstat, [100*alpha/2, 100*(1-alpha/2)]);
    x_pos = x-CI(1,1);
    x_neg = CI(2,1)-x;
    
    y_pos = y-CI(1,2);
    y_neg = CI(2,2)-y;
    
    
    errorbar(x,y,y_neg,y_pos,x_neg,x_pos,'o', ...
        'Color',mc(i, :),"MarkerFaceColor",mc(i, :),'markersize',0.1,'linewidth',2); hold on;
    
    
end
xticks([0:1:20]);
xticklabels(xticks*10);
xlim([5 13])
tch_set_axes;
xlabel('temporal window (ms)');
ylabel('pRF size (deg)');
set(gca,'XMinorTick','on','yMinorTick','on')
title('temporalWindow vs prfSize (CI)')



end