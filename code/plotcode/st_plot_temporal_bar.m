function [fn,saveName] = st_plot_temporal_bar(DT, params,roiList,targetModel,targetVar)
%%
% fn = figure('Renderer', 'painters', 'Position', [10 10 1500 600]);

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
gc = (getmycolors(3));

%%
counter = 1;
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

        for tp = 1:nParam %size(df3.temporal,2)
            agg{sj,eachROI} = df3.(targetVar)(:,tp);
        end
    end
    
end



%%
fn.fig1 = figure('Renderer', 'painters', 'Position', [10 10 1500 600]);
saveName.fig1 = ['subjectwise-mode-median-tau-peak' num2str(targetModel)];

mg = cellfun(@(x) round(x),agg ,'UniformOutput' ,false);
% edges = linspace(0,100,51);
% decc = cellfun(@(x) discretize(x,edges),ecc,'UniformOutput',false);

subplot(121)
plotData = cellfun(@median ,mg);
st_plotBox(plotData,roiList,gc,'\tau (ms)');
title('median');

subplot(122)
plotData = cellfun(@mode ,mg);
st_plotBox(plotData,roiList,gc,'\tau (ms)')
title('mode')


% return;
%%
fn.fig2 = tch_fig(['contributions'],[0,0,40,15]);
saveName.fig2 = ['time-to-peak-streams' num2str(targetModel)];

allKDF=[]; medianDur = [];
gc = (getmycolors(3));
edges = linspace(0,1000);
for er = 1 :length(roiList)
    kdf =[];  kdx=[]; N=[];
    for sj = [1:size(agg,1)]
        data = cell2mat(agg(sj,er));
        if length(data)>1
            [kdf(sj,:) kdx(sj,:)] =ksdensity(data*10, ...
                edges, ...
                'Support','positive','Bandwidth',0.2);
        end
    end
    kdf2= kdf;
    kdf( ~any(kdf,2), : ) = [];  %rows
    [~,mi]=max(kdf');
    medianDur{er} = mi;
    allKDF{er} = kdf2; 
end

nh = hot;
nh(end-25:end,:) = []; % remove white it is too bright
colormap(nh)
% colormap hot

subplot(121);
allRaster =[];
for er = 1:9
    allRaster = [allRaster; normSum(nanmean(allKDF{er}))];
end

imagesc(flipud(allRaster));
xticks([0:5:100])
xticklabels(xticks*10);
yticklabels(fliplr(roiList));
xlabel('Time (ms)'); tch_set_axes; box off; 

xlim([0 30]); % colorbar('XTick', 0:0.25:1);
axis square;   
clip = sort(max(allRaster'));
% caxis([0 mean(clip)+std(clip)]);
caxis([0 0.06]);

colorbar
% caxis([0 clip(end)]);


%%

subplot(122);

gc = (getmycolors(4));
lw = 3;

ev = [allKDF{1}; allKDF{2}; allKDF{3}];
tch_plot_tc(1:length(edges),ev,lw,gc(1,:),gc(1,:),'sem'); hold on;

ev = [allKDF{4}; allKDF{5}];
tch_plot_tc(1:length(edges),ev,lw,gc(2,:),gc(2,:),'sem'); hold on;

ev = [allKDF{6}; allKDF{7}];
tch_plot_tc(1:length(edges),ev,lw,gc(3,:),gc(3,:),'sem'); hold on;

ev = [allKDF{8}; allKDF{9}];

tch_plot_tc(1:length(edges),ev,lw,gc(4,:),gc(4,:),'sem'); hold on;

title('Voxel time-to-peak distribution')
xlim([0 20]);
ylabel('density');
xlabel('time-to-peak(ms)');
axis square;   
xlim([0 35]);
xticks([0:5:30]);
xticklabels(xticks*10)
tch_set_axes;
legend({'Early visual (V1, V2, V3)','Ventral (hV4,VO)','Lateral (LO,TO)','Dorsal (V3AB, IPS)'})
legend boxoff;             
%%



fn.fig3 = figure('Renderer', 'painters', 'Position',...
    [  467         406        1152         343]);
saveName.fig3 = ['ecdf-small-model' num2str(targetModel)];

gc = (getmycolors(3));
for er = length(roiList):-1:1
    v_data = cell2mat(agg(:,er));
    [v_f,v_x]=homemade_ecdf(v_data');
    p = plot(v_x,v_f,'-o','color',gc(er,:),'linewidth',2); hold on;
    p.MarkerFaceColor = gc(er,:);
    p.MarkerSize = 2;
    p.MarkerEdgeColor = gc(er,:);
    xlim([0 50])
    xticks([0:10:100]);
    xticklabels(xticks*10);
    yticks([0:0.25:1])
end
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;


tch_set_axes; box on;
labels = fliplr(roiList);
plots = flipud(get(gca, 'children'));
neworder = [9,8,7,6,5,4,3,2,1];
legend(plots(neworder), labels(neworder),'Fontsize',15)
legend box off;
%%
fn.fig4 = figure('Renderer', 'painters', 'Position',...
    [    463   373   456   373]);
saveName.fig4 = ['ecdf-large-model' num2str(targetModel)];

for er = length(roiList):-1:1
    v_data = cell2mat(agg(:,er));
    [v_f,v_x]=homemade_ecdf(v_data');
    p = plot(v_x,v_f,'-o','color',gc(er,:),'linewidth',2); hold on;
    p.MarkerFaceColor = gc(er,:);
    p.MarkerSize = 2;
    p.MarkerEdgeColor = gc(er,:);
    xticks([0:10:100]);
    xticklabels(xticks*10);
    xlim([1 20]);
    tch_set_axes; box on;

end



%%
fn.fig5 = figure('Renderer', 'painters', 'Position',...
    [       459   142   964   737]);
saveName.fig5 = ['ecdf-indiv-large-model' num2str(targetModel)];

for sj = 1:size(agg,1)
    subplot(5,2,sj)
for er = length(roiList):-1:1
    v_data = cell2mat(agg(sj,er));
    if length(v_data )>2
    [v_f,v_x]=homemade_ecdf(v_data');
    p = plot(v_x,v_f,'-o','color',gc(er,:),'linewidth',1); hold on;
    p.MarkerFaceColor = gc(er,:);
    p.MarkerSize = 2;
    p.MarkerEdgeColor = gc(er,:);
    tch_set_axes; box on;
    xticks([0:10:100]);
    xticklabels(xticks*10);
    xlim([1 30])
    end
end
end


end

function [v_f,v_x] = homemade_ecdf(v_data)
nb_data = numel(v_data);
v_sorted_data = sort(v_data);
v_unique_data = unique(v_data);
nb_unique_data = numel(v_unique_data);
v_data_ecdf = zeros(1,nb_unique_data);
for index = 1:nb_unique_data
    current_data = v_unique_data(index);
    v_data_ecdf(index) = sum(v_sorted_data <= current_data)/nb_data;
end
v_x = [v_unique_data(1) v_unique_data];
v_f = [0 v_data_ecdf];
end
function st_plotBox(plotData,roiList,gc,plotylabel)

b= boxplot(plotData,'Color','k','Notch','off','Widths',0.8);
h = findobj(gcf,'tag','Outliers');
delete(h);
h = findobj(gcf,'tag','Median');
set(h,'color','k')

gc = flipud(gc);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),gc(j,:),'FaceAlpha',.5); hold on;
end
lc = lines(10);
data = (plotData);
n = (length(roiList));
for i = 1:n
    s = scatter(ones(size(data(:,i)))* i ,data(:,i),[], ...
        'MarkerEdgeColor','k','MarkerFaceColor','k', ...
        'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); hold on;
    
%     s = scatter(ones(size(data(:,i)))* i ,data(:,i),[], ...
%         'MarkerEdgeColor',lc(i,:),'MarkerFaceColor',lc(i,:), ...
%         'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); hold on;

end
% yticks = 0:10:100;
yticklabels(yticks*10);
% xticks(mean(reshape(find(nansum(data)),2,[])));
xticklabels(roiList)
h=gca; h. XAxis. TickLength = [0 0];
ylabel(plotylabel);
tch_set_axes;
end



function drawPolar(x0,y0,sigma1,fieldRange,nVoxels,pcolors)
[ph, ecc] = cart2pol(x0, y0);

subX = single(ecc .* cos(ph));
subY = single(ecc .* sin(ph));
subTheta = zeros(size(subX));


% add polar grid on top
p =[];
p.ringTicks = (1:3)/3*fieldRange;
p.gridLineWidth = 2;
p.fontSize = 10;
% p.gridColor = 'k';
polarPlot([], p);

target = 1:length(subX);
% target = 1:nVoxels;

target = Shuffle(target);

if nVoxels < length(subX)
    target = target(1:nVoxels);
end
% finish plotting it
counter = 1;
for i= target
    Center = [subX(i),subY(i)];
    R = sigma1(i);
    h = mycircle(Center,R,1000,'-',pcolors(i,:)); hold on
    h = plot(Center(1),Center(2),'+','MarkerSize',3); hold on
    h.Color =pcolors(i,:);
    counter = counter + 1;
end
end

%% extra ecc and tau
% % % fn.fig0 = figure('Renderer', 'painters', 'Position', [10 10 1500 600]);
% % % edges = linspace(1,12,13);
% % % decc = cellfun(@(x) discretize(x,edges),ecc,'UniformOutput',false);
% % % 
% % % for sj = 1:10
% % % for ttt =1:length(edges)-1
% % %     subplot(1,10,sj)
% % % idx = find(decc{sj,1}==ttt);
% % % pd = agg{sj,1};
% % % if ~isempty(idx)
% % % tch_plot_bar(pd(idx),'r',ttt); hold on;
% % % % xlim([0 10])
% % % end
% % % % tch_plot_bar(pd(idx),[],2); hold on;
% % % end
% % % yticklabels(yticks*10)
% % % ylabel('tau')
% % % xlabel('ecc')
% % % 
% % % end
% % % % agg(sj,er)))
% % % % bar(agg{1,1},ecc{1,1})
% % % saveName.fig0 = ['ecc-tau-' num2str(targetModel) '.png'];

% %%
% fn.fig5 = figure('Renderer', 'painters', 'Position', [     93          31        1796         903]);
% 
% subplot(141)
% lw = 2;
% for er = 1:3
% tch_plot_tc(1:length(edges),allKDF{er},lw,gc(er,:),gc(er,:),'sem'); hold on;
% end
% 
% subplot(142)
% for er = 4:5
% tch_plot_tc(1:length(edges),allKDF{er},lw,gc(er,:),gc(er,:),'sem'); hold on;
% end
% 
% subplot(143)
% for er = 6:7
% tch_plot_tc(1:length(edges),allKDF{er},lw,gc(er,:),gc(er,:),'sem'); hold on;
% end
% 
% 
% subplot(144)
% for er = 8:9
% tch_plot_tc(1:length(edges),allKDF{er},lw,gc(er,:),gc(er,:),'sem'); hold on;
% end
