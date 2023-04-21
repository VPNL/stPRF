function [fn,saveName] = st_compareMethods(DT, params,roiList,targetModel,targetVar)
%%

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
fn.fig1 = figure(1); clf;
saveName.fig1 = ['fmri_t2p_dist'];

% figure(1);
ScrSz = get(0, 'screensize');
set(gcf, 'position',  [ScrSz(1) ScrSz(2) ScrSz(3)/3 ScrSz(4)/2]);

mc = getmycolors(3);
grayyy = [0.25 0.25 0.25];

counter = 1;
for i = [length(roiList):-1:1]
    
    ax = subplot_tight(length(roiList),1, counter,[0.02 0.01]);
    data = (cell2mat(agg(:,i))); hold on
    m = nanmedian(data);
    fmriMedian{i} = m;
    fprintf('median %s: %.2f \n',roiList{i},m*10)
    
    [fi, xi]=ksdensity(data,'support','positive');
    a = area([xi],[fi],'facecolor',mc(i, :),'linewidth',1.5,'edgecolor',grayyy,'FaceAlpha',0.3); hold on;
    plot(m,0,'o','MarkerSize',10,'MarkerEdgecolor',mc(i, :),'MarkerFacecolor','w','linewidth',4); hold on
    line([0 5],[0 0],'color',grayyy,'linewidth',1.5); hold on;
    
    if i > 1
        set(ax,'XLim',[0,30],'ylabel',[],'ytick',[],'ycolor','none', ...
            'xlabel',[],'xtick',[],'xcolor','w', 'tickDir','out'); box off;
        
    else
        set(ax,'XLim',[0,30],'ylabel',[],'ycolor','none', ...
            'xaxisLocation','origin','xticklabel',[],'linewidth',2, ...
            'XMinorTick','on','TickDir','out','xcolor','none'); box off;
        
    end
    counter = counter+1;
end

% printnice(gcf,[1 300],plotDir,['fmri_t2p_dist']);
% 


%% greon 2022
g = load('./results/paperfMRI/other/groen.mat');

fn.fig2 = figure(2); clf;
saveName.fig2 = ['ecog_groen_t2p'];
ScrSz = get(0, 'screensize');
set(gcf, 'position',  [ScrSz(1) ScrSz(2) ScrSz(3)/3 ScrSz(4)/2]);

mc = getmycolors(3);
grayyy = [0.25 0.25 0.25];

counter = 1;
for er = [length(roiList):-1:1]
    %     er
    ax = subplot_tight(length(roiList),1, counter,[0.02 0.01]);
    
    [idx] = ismember(upper(g.areaNames),roiList{er});
    
    if isempty(find(idx, 1))
        set(ax,'XLim',[0,250] ,'ytick',[],'ycolor','none', ...
            'xtick',[],'tickDir','out','Fontsize',14, ...
            'linewidth',2,'xcolor','none'); box off;
        counter = counter +1;
        continue;
    end
    
    
    data = cell2mat(g.t2p(:,idx))*1000;
    m = nanmedian(data);
    [fi, xi]=ksdensity(data,'support','positive');
    a = area([xi],[fi],'facecolor',mc(er, :),'linewidth',1.5,'edgecolor',grayyy,'FaceAlpha',0.3); hold on;
    plot(m,0,'diamond','MarkerSize',12,'MarkerEdgecolor','w','MarkerFacecolor',mc(er, :),'linewidth',2); hold on
    line([0 min(xi)],[0 0],'color','k'); hold on;
    line([max(xi) 1000],[0 0],'color','k');
    
    plot(fmriMedian{er}*10,0,'o','MarkerSize',10,'MarkerEdgecolor',mc(er, :),'MarkerFacecolor','w','linewidth',4); hold on
    set(ax,'XLim',[0,300] ,'ytick',[],'ycolor','none', ...
        'xtick',[],'tickDir','out','Fontsize',14, ...
        'linewidth',2,'xcolor','none'); box off;
    counter = counter+1;

end



%%

fn.fig3 = figure(3); clf;
saveName.fig3 = ['electro_t2p'];
ScrSz = get(0, 'screensize');
set(gcf, 'position',  [ScrSz(1) ScrSz(2) ScrSz(3)/3 ScrSz(4)/2]);

e = load('./results/paperfMRI/other/electro.mat');
T = e.T;

counter = 1;
mROI = {'V1','V2','V3','V4','TEm/TEa','','MT','V3A' 'LIP'};
locs = [1,     2,  3,   4    5, 6,  7,  8 , 9 ];
yt = cell(1,length(mROI));
mc=getmycolors(3);
for er = 1 :length(mROI)
    idx = T.target_area == mROI{er};

    if sum(idx) > 0
        t = T(idx,:);
        if size(t.latency,1) > 1
            electro = nanmean(t.latency);
        else
            electro = (t.latency);
        end
    else
        electro = [NaN NaN NaN];
    end

    errorbar(electro(2),locs(er), 0,0,electro(2)-electro(1),electro(3)-electro(2),'diamond',...
        "MarkerFaceColor",mc(locs(er), :),'Color',mc(locs(er), :),'markersize',10,'linewidth',2 ); hold on
    
    v_data = median(cell2mat(agg(:,er))); hold on
    set(gca,'XLim',[0,300],'yLim',[0 10],'ycolor','w','ytick',[],'ylabel',[], ...
        'XMinorTick','on');
    
    plot(fmriMedian{er}*10,locs(er),'o','MarkerSize',8,'MarkerEdgecolor',mc(er, :),'MarkerFacecolor','w','linewidth',2.5); hold on
    counter = counter+1;
    yt(locs(er)) = mROI(er);
end
yticks(locs)
yticklabels(mROI)
tch_set_axes;

%% extra

%% Zhou 2019
% % % % fn.fig4 = figure(4); clf;
% % % % saveName.fig4 = ['ecog_zhou_t2p'];
% % % % 
% % % % b = load('./results/paperfMRI/other/dn_params.mat');
% % % % dn   = b.prm.ecog.dn;
% % % % t2pk(1, :) = dn.summaryMetrics.t2pk; % time to peak
% % % % 
% % % % nRois = 3;
% % % % pt= [10 50 90];
% % % % counter = 1;
% % % % for er = nRois:-1: 1
% % % %     ax = subplot_tight(length(roiList)+1,1, 11-er,[0.015 0.02]);
% % % % 
% % % %     idx = (er - 1) * 100 + 1 : er * 100;
% % % %     v_data    = t2pk(1, idx);
% % % %     s =  prctile(v_data,pt);
% % % %     [kdf kdx] =ksdensity(v_data);
% % % %     
% % % %     a = area(kdx,kdf,'facecolor',mc(er, :),'linewidth',2,'edgecolor',[0.76 0.75 .75]); hold on
% % % %     a.FaceAlpha = 0.2;
% % % %  
% % % %     plot(s(2),0,'o','MarkerSize',10,'MarkerEdgecolor','w','MarkerFacecolor',mc(er, :)); hold on
% % % %     xlim([0 400])
% % % %     set(ax,'XLim',[0,300],'ylabel',[],'ytick',[],'ycolor','w','xlabel',[],'xtick',[],'xcolor','w'); box off;
% % % %     counter = counter+1;
% % % % % % % % % % end
% % % % % % 
% % % % 
%% fMRI BOOT STRAP
% % % % 
% % % % fn.fig2 = figure(2); clf;
% % % % saveName.fig2 = ['fmri_t2p_boot'];
% % % % 
% % % % 
% % % % ScrSz = get(0, 'screensize');
% % % % set(gcf, 'position',  [ScrSz(1) ScrSz(2) ScrSz(3)/3 ScrSz(4)/2]);
% % % % mc = getmycolors(3);
% % % % grayyy = [0.25 0.25 0.25];
% % % % rng('default'); clc; counter = 1;
% % % % for er = [length(roiList):-1:1]
% % % %     
% % % %         ax = subplot_tight(length(roiList),1, counter,[0.02 0.01]);
% % % %         roi_data = agg(:,er);
% % % %         fun = @(x) repmat(find(cellfun(@(y) isequal(x,y), roi_data)), 1, length(x));
% % % %         C_new = cellfun(fun, roi_data, 'UniformOutput', false);
% % % %         C_new = C_new(~cellfun('isempty',C_new));
% % % %         strata = cell2mat(cellfun(@(x) permute(x, [2 1]), C_new, 'UniformOutput', false));
% % % % 
% % % % 
% % % %         data = cell2mat(roi_data);
% % % % 
% % % %         % get the unique elements and compute the proportion of each unique element
% % % %         [~,~,unique_indices] = unique(strata);
% % % %         proportions = accumarray(unique_indices,1)./length(strata);
% % % %         % replace each element in the original array with its computed proportion
% % % %         wts = arrayfun(@(x) proportions(unique_indices(x)), 1:length(data));
% % % % 
% % % % %         bootstat = bootstrp(100,@median,data,'weights',wts);
% % % %         bootstat = bootstrp(100,@median,data);
% % % % 
% % % %         ci = prctile(bootstat,[25.5 50  97.5]);
% % % %         [fi, xi]=ksdensity(bootstat);
% % % %         
% % % %         % plot
% % % %         a = area([xi],[fi],'facecolor',mc(er, :),'linewidth',2,'edgecolor',grayyy,'FaceAlpha',1); hold on;
% % % %         scatter(mean(bootstat),0,50,'MarkerEdgeColor','w','MarkerFaceColor',mc(er, :),'Linewidth',1); hold on;
% % % % 
% % % %         line([0 min(xi)],[0 0],'color',grayyy); hold on;
% % % %         line([max(xi) 100],[0 0],'color',grayyy);
% % % %         
% % % %         if  er  == 1
% % % %             set(ax,'XLim',[0,30] ,'ytick',[],'ycolor','none', ...
% % % %                 'tickDir','out','XMinorTick','on','Fontsize',14, ...
% % % %                 'linewidth',2,'xcolor','k'); box off; 
% % % %         else
% % % %             set(ax,'XLim',[0,30] ,'ytick',[],'ycolor','none', ...
% % % %             'xtick',[],'tickDir','out','Fontsize',14, ...
% % % %             'linewidth',2,'xcolor','none'); box off; 
% % % %         end
% % % % 
% % % % 
% % % %         
% % % %         counter = counter + 1;
% % % % 
% % % %     
% % % % end

%%
% 
% 
% paper = {'Schmolesky';'Nowak'; 'Celebraini';'Maunsell'; ...
%     'Knierim';'Vogels';'Celebrini'; 'Raiguel';'Bair'};
% year = [1998;1995;1993;1992; ...
%     1992;1994;1993;1989;2002];
% target_area = ["V1";"V1";"V1";"V1";...
%         "V1";"V1";"V1";"V1";"V1"];
% latency = [51 66 77;   45 77 187; 41 62 151;  27 45 65; ...
%            42.5 55 75; 40 60 120; 47.5 60 90; 45 85 267; 26 52 87];
% T1 = table(paper,year,target_area,latency);
% 
% paper = {'Schmolesky';'Nowak'; 'Raiguel';'Schmolesky';'Schmolesky'...
%     ;'Zamarashkina';'Chang'};
% year = [1998; 1995;1989;1998;1998; ...
%     2020;2014];
% target_area = ["V2";"V2";"V2";"V3";"V4";...
%     "V4";"V4"];
% latency = [64 82 104; 56 84.5 120; 71 96 185;60 72 81; 74 104 159; ...
%     25 75.5 204; 60 85 150];
% T2 = table(paper,year,target_area,latency);
% 
% 
% paper = {'Schmolesky';'Raiguel';'Raiguel';'Bair';'Nakhla';'Nakhla'};
% year = [1998;1989;1999;2002;2021;2021];
% target_area = ["MT";"MT";"MT";"MT";"MT";"V3A"];
% latency = [58 72 85;35 94 273;35 87 325 ;29 40 48;40 74 195; 40 78 243];
% T3 = table(paper,year,target_area,latency);
% 
% 
% 
% paper = {'Schmolesky';'Thompson';'Bushnell';'Baylis';'Perrett';'Barash'};
% year = [1998;1996;1991;1987;1982;1991];
% target_area = ["FEF";"FEF";"FEF";"TEm/TEa";"STS";"LIP"];
% latency = [54 75 91; 48 67 106; 40 67 118; 89 112 176; 90 110 143;70 100 180];
% T4 = table(paper,year,target_area,latency);
% 
% 
% T = [T1 ;T2; T3; T4];
% 

end

