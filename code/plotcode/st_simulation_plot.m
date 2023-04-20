function st_simulation_plot(gt,pt)

% gt: ground-truth table
% pt: predicted table (model estimated)


% simDir    = '/share/kalanit/users/insubkim/oak/biac2/kgs/projects/spatiotemporal/experiments/simulator';
% VoxelName = 'voxel-all_noise2-3i50f';

% effectivesigma = 0;
% 

% cd(simDir)
% sessionDir = fullfile(simDir,'data','subj01');
% noiselevel = split(VoxelName,'_');
% noiselevel = noiselevel{2};
% noiselevel = split(noiselevel,'-');
% noiselevel = noiselevel{1};

%   '081422'
% analysisDate = Constants.getDir.sessionDate;
% pt = load(fullfile(sessionDir,'dataTable', ...
%     Constants.getDir.sessionDate, ...
%     'gray',VoxelName,'df_subj01.mat'));


tmp = [];
tmp.tmodel   = pt.DT.tmodel;
tmp.testData = cell2mat(pt.DT.testData);
tmp.testPred = cell2mat(pt.DT.testPred);
tmp = struct2table(tmp);

params = [];
params.tmodel =  "spatial";
tc{1} = tableFilter(tmp,params);
params.tmodel =  "DN_ST";
tc{2} = tableFilter(tmp,params);
params.tmodel = 'CST';
tc{3} = tableFilter(tmp,params);
pt.DT= removevars(pt.DT,{'trainData','testData','trainPred','testPred','trainSet','testSet'});


% gt = getAllFiles(fullfile('./results/',analysisDate,'voxels/'),...
%     sprintf('GT*%s.mat',noiselevel),1);
% gt = load(gt{1});
gt = gt.DT;


%fix order of the ground-truth data
current_order = {gt{1}.tmodel{1}, gt{2}.tmodel{1}, gt{3}.tmodel{1}};

% Define the desired order
desired_order = {'spatial', 'DN-ST', 'CST'};

% Find the indices of each element in the original array
[~, idx] = ismember(desired_order, current_order);

gt = gt(idx);
targetTmodel =  desired_order;


DT = pt.DT;

% give ECC values
DT.ecc = pdist2([0 0], [DT.x0,DT.y0])';
rows = find(ismember(DT.tmodel, {'1ch_dcts', 'ST_DN'}));

% rename your models
DT.modelNumber( find(ismember(DT.tmodel, {'1ch_glm', 'spatial'})))  = 1;
DT.modelNumber( find(ismember(DT.tmodel, {'1ch_dcts', 'DN_ST'}))) = 2;
DT.modelNumber(  find(ismember(DT.tmodel, {'3ch_stLN', 'CST'}))) = 3;

% effectivesigma = 0;
% if effectivesigma
%     DT((DT.tmodel=="CST"),:).sigma = ...
%         DT((DT.tmodel=="CST"),:).sigma./sqrt(DT((DT.tmodel=="CST"),:).exponent);
%     gt{3}.RF(:,3)
%     gt3_exp = cell2mat(gt{3}.tparam);
%     gt3_exp = gt3_exp(:,3);
%     gt{3}.RF(:,3) = gt{3}.RF(:,3)./sqrt(gt3_exp);
% 
% end


params = []; df=[]; df_all=[]; pt=[];
for tm = 1:length(targetTmodel)
    params.modelNumber = tm;
    df = tableFilter(DT,params);
    df = tableArrayFormat(df);
    
    if length(unique(DT.split)) > 1
        df= varfun(@mean,df,'InputVariables',...
            {'x0','y0','sigma','exponent','temporal','cv_varexp','beta','varexp','searchFit'} ,...
            'GroupingVariables',{'coords','modelNumber'});
        df = renamevars(df,df.Properties.VariableNames, ...
            strrep(df.Properties.VariableNames,'mean_',''));
    end
    df = mergevars(df,{'x0','y0','sigma'},'NewVariableName','RF');
    
    
    df = renamevars(df,'temporal', ...
        'tparam');
    
    pt{tm} = df;
    
end

tmn  = {'spatial','DN-ST','CST'};
ts   = {'\itx','\ity','\sigma'};
tms1 = {'\tau_1', '\tau_2', '\itn_D_N', '\sigma_D_N'};
tms2 = {'\tau', '\itn'};


% 3 6 35 54 91 72
%% SNR and Varexp
mycolor = getmycolors(1);

% fn = figure('units','centimeters','position',[0,0,10,10]);
tch_fig('SNR-VAREXP',[0,0,20,10]);
for tm = 1:length(targetTmodel)
    subplot(121)
    g = gt{tm};
    p = pt{tm};
    xgrid = linspace(-20,20,1000)';
    if sum(isinf(g.SNR) ) < 0
        pd = fitdist(g.SNR,'Normal');
        pdfEst = pdf(pd,xgrid);
        line(xgrid,pdfEst,'color',mycolor(tm,:),'linewidth',3);
    %     xlim([-0.25 0.4]); axis square;
        axis square;
        ylabel('density');
        xlabel('SNR (dB)');
        tch_set_axes;
    end
    
    subplot(122);
    binWidth = 0.01;
    xgrid = linspace(0,1,1000)';
    lastVal = ceil(max(p.cv_varexp));
    pd = fitdist(p.cv_varexp,'Normal');
    pdfEst = pdf(pd,xgrid);
    line(xgrid,pdfEst,'color',mycolor(tm,:),'linewidth',3); hold on;
%     xlim([0.25 0.75]); axis square;
%     xticks(0:0.25:1)
    xlabel('model accuracy (R^2)');
    ylabel('density');
    tch_set_axes;
    


end
% legend(tmn,'location','best'); legend box off;

%% plot scatter
tch_fig('xys',[0,0,40,20]);
counter = 1;
Fsize = 10;
pMargin = [0.01,0.05];
pMargin = [0.01,0.03];

for tm = 1:3
    g = gt{tm};
    p = pt{tm};
        
    gtmp = cell2mat(g.tparam);
    if tm == 1
        tn = [ts];
        gtmp = [];
        ptmp = [];
    elseif tm == 2
        tn = [ts tms1];
        gtmp = gtmp(:,[1,3,4,5]);
        ptmp = p.tparam(:,[1,3,4,5]);
    elseif tm ==3
        tn = [ts tms2];

       gtmp = gtmp(:,[1,3]);
       gtmp(:,1) = gtmp(:,1)./100;
       ptmp = [p.tparam(:,1)./100 p.exponent];

    end
    
    g.params = [g.RF gtmp];
    p.params = [p.RF ptmp];

    % calculate errors
    mSTP =[]; stdSTP=[];
    eSTP{tm} = abs(g.params - p.params);
    mSTP(tm,:) = round(mean(eSTP{tm}),2);
    stdSTP(tm,:) = round(std(eSTP{tm}),2);

    
    for s = 1:size(g.params,2)
       if tm ==2 && s == 1
          counter = 8;      
       end
           
         
        subplot_tight(3,7,counter,pMargin);
        scatter(g.params(:,s),p.params(:,s),'MarkerEdgeColor',mycolor(tm,:)); hold on;
        f=fit(double(g.params(:,s)),double(p.params(:,s)),'poly1');



            if s ==1
                ylabel([tmn{tm},' prediction']);
            else
                ylabel({''});
                
            end
        
        if s <3
            xy = [-10 10];
            xlim([-10 10])
            ylim([-10 10])
            xticks([-10:5:10]);
            yticks([-10:5:10]);
        elseif s == 3
            xy = [0 5];
            xlim([0 4])
            ylim([0 4])
            xticks([0:4]);
            yticks([0:4]);
        else
            if tm ==2
                if s == 4 || s == 5
                    TickL = [0 0.5 1];
                elseif  s == 6
                    TickL = [1 3.5 6];
                elseif s ==7 
                    TickL = [0 0.25 0.5];
                end
                xticks(TickL);
                yticks(TickL);
                xlim([min(TickL) max(TickL)]);
                ylim([min(TickL) max(TickL)]);
                xy = [min(TickL)  max(TickL)];

            elseif tm == 3
                TickL = [0 0.5 1];
                xticks(TickL);
                yticks(TickL);
                xlim([min(TickL) max(TickL)]);
                ylim([min(TickL) max(TickL)]);
                xy = [min(TickL)  max(TickL)];
            end
            
        end

        line(xy,xy,'Color','k','LineWidth',1.5,'LineStyle','--'); hold on;
        title(sprintf('%s',tn{s})); 
        axis square; box off; tch_set_axes; legend off;

        counter = counter + 1;
    end
end
% printnice(5,[1 300],saveDir,sprintf('scatter-params'));

%%  absoulte percentage error (APE)
mycolor = getmycolors(1);

% tch_fig('MAPE',[0,0,36, 9]);
tch_fig('APE',[0,0,20, 5]);
thresh = 90; 

ae = @(y,yhat) (abs(y-yhat));
ape = @(y,yhat) (abs((y-yhat)./y))*100;

ep=[]; ee=[];
for tm = 1:3
    g = gt{tm};
    p = pt{tm};
    gtmp = cell2mat(g.tparam);
    if tm == 1
        gtmp = [];
        ptmp = [];
    elseif tm == 2
        tn = [ts tms1];
        gtmp = gtmp(:,[1,3,4,5]);
        ptmp = p.tparam(:,[1,3,4,5]);
    elseif tm ==3
        tn = [ts tms2];
        gtmp = gtmp(:,[1,3]);
        gtmp(:,1) = gtmp(:,1)./100;
        ptmp = [p.tparam(:,1)./100 p.exponent];
    end
    g.params = [g.RF gtmp];
    p.params = [p.RF ptmp];
    
    for s = 1:size(g.params,2)
        y = g.params(:,s);
        yhat = p.params(:,s);
        
        [R,P] = corr(y,yhat,'Type','Pearson') ;
        
        cr{tm,s} = R;
        cp{tm,s} = P;
        y(y == 0) = 0.0001;
        ep{tm,s} = ape(y,yhat);
        ee{tm,s} = ae(y,yhat);
        
        
    end
end
clc;
disp("correlation")
disp(cellfun(@(x) round(x,3), cr, 'UniformOutput', false))


TF = cellfun(@(x) ~isoutlier(x, 'percentiles', [0 thresh]), ep, 'UniformOutput', false);
ep = cellfun(@(x,y) x(y), ep, TF, 'UniformOutput', false);


disp("median percent error")
disp(cellfun(@(x) round(median(x),3), ep, 'UniformOutput', false))


plotData =ep(1:3,1:3);
plotData = plotData(:);
plotData = cell2mat(plotData');

plotData_gap = zeros(size(plotData,1),11);
plotData_gap(:,[1:3 5:7 9:11]) = plotData;
plotData = plotData_gap;

% DN params
plotData =ep(2,4:end);
plotData = plotData(:);
plotData = cell2mat(plotData');

%     plotData = removeOutlier(plotData,thresh);
plotData_gap = [plotData_gap  zeros(size(plotData,1),1) plotData];

% CST params

plotData =ep(3,4:5);
plotData = plotData(:);
plotData = cell2mat(plotData');

%     plotData = removeOutlier(plotData,thresh);
plotData_gap = [plotData_gap  zeros(size(plotData,1),1) plotData];

g_color = repmat([mycolor; 1 1 1 ],3,1);
g_color = [g_color; repmat(mycolor(2,:),4,1);[0 0 0];repmat(mycolor(3,:),2,1)];


groupm = zeros(size(plotData_gap));
for n= 1:size(plotData_gap,2)
    groupm(:,n) =    ones(size(groupm,1),1)*n;
end

vs = violinplot(plotData_gap, ...
    groupm, ...
    'ViolinColor',g_color, ...
    'ViolinAlpha', .8, ...
    'MedianColor',[0 0 0], ...
    'BoxColor',[.9 .9 .9], ...
    'ShowData',false, ...
    'ShowBox',false, ...
    'QuartileStyle','boxplot', ...
    'ShowWhiskers',false);


for n= [4,8,12,17]
    vs(1,n).ShowData = 0;
    vs(1,n).ShowMedian = 0;
    
end

% xlabel({'First line';'Second line'})

%     tch_plot_bar(plotData_gap,g_color);
ylabel({'absolute';'percentage error'})
xticks([2 6 10  13 14 15 16 18 19]); yticks([0:25:100])
xticklabels([ts tms1 tms2]); tch_set_axes;
ax = gca;

ax.YAxis.FontSize = 12; %for y-axis ;
ax.XAxis.FontSize = 10; %for y-axis ;
set(gca,'TickDir','in');
ylim([0 110])
%     h = legend(tmn,'location','northeast'); legend box off
% printnice(ff,[1 300],saveDir,sprintf('fig-%d',ff));
% printnice(ff,0,saveDir,sprintf('fig-%d',ff));

%% tc plot
s=1;
for tm = 1:3
    g = gt{tm};
    pp = [g.RF(s,:) g.tparam{1,:}];
    
    if tm ==1
        pp = pp(1:3);
        tl = sprintf('%s: %0.2f %s: %0.2f %s: %0.2f\n    ', ...
            ts{1},pp(1),ts{2},pp(2),ts{3},pp(3));
        tl = sprintf('LS:[%0.2f, %0.2f, %0.2f]', ...
            pp(1),pp(2),pp(3));
        
        
    elseif tm == 2
        pp(5) =[];
        %         tl = sprintf('%s: %0.2f %s: %0.2f %s: %0.2f\n%s: %0.2f %s: %0.2f %s: %0.2f %s: %0.2f %s: %0.2f', ...
        %             ts{1},pp(1),ts{2},pp(2),ts{3},pp(3), ...
        %             tms1{1},pp(4), tms1{2},pp(5),tms1{3},pp(3),tms1{4},pp(7));
        tl = sprintf('DN-ST:[%0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f]', ...
            pp(1),pp(2),pp(3),pp(4),pp(5),pp(6),pp(7));
        
    elseif tm == 3
        pp(4) =[];
        %        tl = sprintf('%s: %0.2f %s: %0.2f %s: %0.2f\n%s: %0.2f %s: %0.2f', ...
        %            ts{1},pp(1),ts{2},pp(2),ts{3},pp(3), ...
        %            tms2{1},pp(4), tms2{2},pp(5));
        tl = sprintf('CST:[%0.2f, %0.2f, %0.2f, %0.2f, %0.2f]', ...
            pp(1),pp(2),pp(3),pp(4)./100,pp(5));
        
    end
    tls{tm} = tl;

end

tch_fig('tc-plot-simple',[0,0,35,15]);
counter = 1;

for tm = 1:3
    g = gt{tm};
    ptc = tc{tm};

    % 5 15 52 53  88
    % 73
    for nl = 1:2
        subplot(2,1,nl)
        
        if nl == 1
            plot((g.tc(s,:)), ...
                'linewidth',3,'color',mycolor(tm,:)); hold on;

        elseif nl ==2
            plot((ptc.testData(s,:)), ...
                'linewidth',3,'color',mycolor(tm,:)); hold on;
            xlabel('time (s)', 'fontsize', 10);
        end
        ylabel('% signal change', 'fontsize', 10);
        xlim([ 0   420]); 
        ylim([-7 15])
        tch_set_axes;
        counter = counter + 1;

    
        
    if nl ==1
        if tm == 1
            ht = 19;
        elseif tm == 2
            ht = 16.6;
        elseif tm == 3
            ht = 14;

        end
        text(250,ht,tls{tm},'FontSize',10,'BackgroundColor','w', ...
        'HorizontalAlignment','left','Color',mycolor(tm,:), ...
        'FontSize',10,'FontWeight','Bold'); hold on
    end
%         dim = [.7 .8 .2 .2];
%         h = annotation('textbox', dim, 'string', tls, ...
%             'FitBoxToText','on','Color','w')

    end

end
end

