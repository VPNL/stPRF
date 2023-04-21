%%
function [fn,saveName]=st_plot_beta(DT, params,roiList)
% params=[];
% initial filter
df = tableFilter(DT,params);

% % params=[]; df=[];
% targetTmodel =  {'1ch_glm','1ch_dcts', '3ch_stLN'};
% for tm = 1:3
%     params.tmodel = targetTmodel{tm};
%     df{tm} = tableFilter(DT,params);
%     df{tm} = tableOutlier(df{tm},'beta');
% end
% df = [df{1}; df{2}; df{3}];

%%



% give number to ROIs
df = st_tablegiveROInumber(df, roiList);

params=[];
params.modelNumber = 3;
df2 = tableFilter(df,params);

% % unpack table cells to arrays
df2 = tableArrayFormat(df2);

ecc = 0:12;

% add Sus/trans logical
sus =  df2.beta(:,1) > df2.beta(:,2) ;
ratio = df2.beta(:,1) ./ df2.beta(:,2) ;
eccBin = discretize(df2.ecc,[0:2:12]);
nBin = length(unique(eccBin));

df2 = addvars(df2,sus);
df2 = addvars(df2,ratio);
df2 = addvars(df2,eccBin);

%%
bv = [];
for eachROI =1:length(roiList)
    params.roiNumber = eachROI;
    df3= tableFilter(df2,params);
    G = groupsummary(df3,{'subjID'}, ...
        'mean','beta','IncludeEmptyGroups',true);
    if ~isempty(setdiff(unique(df2.subjID), G.subjID))
        newG = []; % counter = 1;
        ms = (setdiff(unique(df2.subjID), G.subjID));
        for i = 1:length(ms)
            newG.subjID = ms(i);
            newG.GroupCount = 0;
            newG.mean_beta  = [NaN NaN];
            G = [G ; struct2table(newG)];

        end
    end
    bv{eachROI} = G.mean_beta;
end


%%
fn.fig1 = tch_fig(['contributions'],[0,0,10,10]);
saveName.fig1 = ['contributions_oct'];

x = linspace(0,3,10);
y = x;
plot(x,y,'k--')

mc = getmycolors(3);
for eachROI =1:length(roiList)

betas = double(cell2mat(bv(eachROI)));
betas = betas(~isnan(betas(:,1)),:);

means = nanmean(betas);
sems = nanstd(betas) / sqrt(size(betas, 1) - 1);
ncats = 1;
fill_cols = [.5 .5 1; 1 .5 .5; .5 .5 .5];
line_cols = [0 0 1; 1 0 0; 0 0 0];
xm = means(:, ncats + 1:ncats * 2);
ym = means(:, 1:ncats);
xe = sems(:, ncats + 1:ncats * 2);
ye = sems(:, 1:ncats);

xmax = max(xm + xe); ymax = max(ym + ye);
xmin = min(xm - xe); ymin = min(ym - ye);
lims = [min([0 xmin ymin]) max([ceil(xmax) ceil(ymax)])];


oo = draw_oct(xm, ym, xe, ye, mc(eachROI, :), 'k'); hold on;
text(xm , ym, roiList{eachROI}, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'FontName', 'Helvetica','Color','w','FontWeight','BOLD');
end

leftshift = 0.06;
draw_oct(0.55, 1.9, 0.05, 0.05, [0 0 1], 'k'); hold on;
text(0.55+leftshift , 1.9, 'Early visual', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'FontName', 'Helvetica');

draw_oct(0.55, 1.8, 0.05, 0.05, [0.1 0.5 0.3], 'k'); hold on;
text(0.55+leftshift , 1.8, 'Ventral', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'FontName', 'Helvetica');

draw_oct(0.55, 1.7, 0.05, 0.05, [1 0.5 0], 'k'); hold on;
text(0.55+leftshift , 1.7, 'Lateral', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'FontName', 'Helvetica');

oo = draw_oct(0.55, 1.6, 0.05, 0.05, [1 0 1], 'k'); hold on;
text(0.55+leftshift , 1.6, 'Dorsal', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'FontName', 'Helvetica');

axis square; xticks([0.5 1 1.5 2]); yticks([0.5 1 1.5 2])
tch_set_axes; xlim([0.47 3]); ylim([0.47 3]);
xlabel('Transient \beta_T');
ylabel('Sustained \beta_S');

%%
allBs = []; allBt = [];
SN = (unique(df2.subjID));
for es = 1:length(SN)
    for eachROI =1:length(roiList)
        params.roiNumber = eachROI;
        params.subjID    = SN(es);
        df3= tableFilter(df2,params);
        y_s =[]; y_t =[];
        if length(df3.ecc) > 30
            nsample = 30;
            nBoot = 100;
        else
            nsample =length(df3.ecc);
            nBoot = 0;
            y_s = NaN(size(ecc));
            y_t = NaN(size(ecc));
        end
        for eachboot = 1:nBoot
            i = randsample(length(df3.ecc),nsample);
            [y_s(eachboot,:)]= p2fit(df3.ecc(i),df3.beta(i,1),ecc);
            [y_t(eachboot,:)]= p2fit(df3.ecc(i),df3.beta(i,2),ecc);
        end
        y_s = nanmean(y_s);
        y_t = nanmean(y_t);
        allBs{es,eachROI} =y_s;
        allBt{es,eachROI} =y_t;
    end
end
%%
fn.fig2 = tch_fig(['ecc_beta'],[0,0,10,14]);
saveName.fig2 = ['ecc_beta'];

plot_locations = [ 1 2 3 4 5 7 8 10 11];

    for eachROI =1:length(roiList)
        idx = find(cell2mat(cellfun(@(x) sum(isnan(x)), allBs(:,eachROI), 'UniformOutput', false))==0);
        subplot_tight(4,3,plot_locations(eachROI),[0 0.1])
       tch_plot_tc(ecc,cell2mat(allBs(idx,eachROI)),2,[0 0 0],[0 0 0],'sem');
       tch_plot_tc(ecc,cell2mat(allBt(idx,eachROI)),2,[1 0 0],[0 0 0],'sem');

        [~,sMax(eachROI)]= max(mean(cell2mat(allBs(idx,eachROI))));
        [~,tMax(eachROI)]= max(mean(cell2mat(allBt(idx,eachROI))));
        
        title(roiList(eachROI));
        xticks(0:3:12);xticklabels([0:3:12]); axis square;
        ylim([0 2.5]); xlim([0 12]); tch_set_axes; 
    end
   



%%
fn.fig3 = tch_fig(['betaRatio'], [0,0,15,7]);
saveName.fig3 = ['betaRatio_bar'];

plotData = cellfun(@(x) x(:,1)./x(:,2),bv,'UniformOutput',false);

mc = flipud(getmycolors(3));

% plotData = ysustrans;
b= boxplot([cell2mat(plotData)],'Color','k','Notch','off','Widths',0.7);
h = findobj(gcf,'tag','Outliers');
delete(h);
h = findobj(gcf,'tag','Median');
set(h,'color','k')

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),mc(j,:),'FaceAlpha',.5); hold on;
end

data =cell2mat(plotData);
n = (length(roiList));
for i = 1:n
    s = scatter(ones(size(data(:,i)))* i ,data(:,i),[], ...
            'MarkerEdgeColor','k','MarkerFaceColor','k', ...
            'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); hold on;
end
h =[]; p =[];
for er = 1:length(roiList)
    betas = double(cell2mat(plotData(er)));
    betas = betas(~isnan(betas(:,1)),:);
    [h{er}, p{er}]=ttest(betas,1);
end


h = cell2mat(h);
p = cell2mat(p);
for i = 1:n
    if h(i) ==1
        if p(i) < .001
            s = text( i ,2.3,'***','FontSize',12); hold on;
        elseif p(i) < .01
            s = text( i ,2.3,'**','FontSize',12); hold on;
        else
        
            s = text( i ,2.3,'*','FontSize',12); hold on;
        end
    end
end
% 
ylim([0 2.5]); yticks([0:0.5:2.5]); yticklabels([0:0.5:2.5]);
xticklabels(erase(roiList,'rh_'))
h=gca; h. XAxis. TickLength = [0.01/4 0.025]; 
ylabel([{'contribution ratio'},{'(\beta_s / \beta_t)'}]);
tch_set_axes;


end

