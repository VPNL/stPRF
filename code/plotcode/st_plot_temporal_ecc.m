function [fn,saveName] = st_plot_temporal_ecc(DT, params,roiList,targetModel)
%%
params =[];
% fn = figure('Renderer', 'painters', 'Position', [100 100 500 500]);

% modelNumber = targetModel;
params.modelNumber  = targetModel;
% initial filter
df = tableFilter(DT,params);

% give number to ROIs
df = st_tablegiveROInumber(df, roiList);

thresh = [];
thresh.temporal = 99;
df = tableThreshold(df,thresh,'upper'); % remove anything upper;
% 
thresh = [];
thresh.temporal = 4;  % max ecc
df = tableThreshold(df,thresh,'lower'); % remove voxels above 12 ecc

gc = (getmycolors(3));


edges=[0:1:12];
df.eccBin = discretize(df.ecc,edges);

all=[];
params=[];
params.modelNumber = targetModel;
df2 = tableFilter(df,params);
df2 = tableArrayFormat(df2);
id = unique(DT.subjID);
T_eccTemporal =[]; T_eccSigma=[];




%%
binsize = 2;
% % unpack table cells to arrays
df2 = tableArrayFormat(df2);
id = unique(DT.subjID);
T_eccTemporal =[]; 
T_eccSigma=[];
for eachROI = 1:length(roiList)
    allTable=[];
    for sj = 1:length(id)
        params.subjID =  id(sj);
        params.roiNumber = eachROI;
        df3 = tableFilter(df2,params);

        if ~isempty(df3) && height(df3) > 20    
            data{sj,eachROI} =st_plot_EccTemporal(df3,binsize);
            T_eccTemporal = table(repmat(sj,size(data{sj,eachROI}.x)), ...
                data{sj,eachROI}.x, ...
                data{sj,eachROI}.y, ...
                'VariableNames',{'subject','ecc','temporal'});

            data2{sj,eachROI} =st_plot_EccSigma(df3,binsize);
            T_eccSigma = table(repmat(sj,size(data2{sj,eachROI}.x)), ...
                data2{sj,eachROI}.x, ...
                data2{sj,eachROI}.y, ...
                'VariableNames',{'subject','ecc','sigma'});
            eccTable = join(T_eccTemporal,T_eccSigma);
        else
            % empty 
            cdv=3;

        end
        
    allTable = [allTable; eccTable];
    roiTable{eachROI} = allTable;
    end
    
    lme_temporal_null{eachROI} = fitlme(roiTable{eachROI}, 'temporal ~ ecc');
    lme_temporal{eachROI} = fitlme(roiTable{eachROI}, 'temporal ~ ecc  + (1|subject)');
    lme_temporal2{eachROI} = fitlme(roiTable{eachROI}, 'temporal ~ ecc  + (ecc|subject)');

    lme_sigma{eachROI} = fitlme(roiTable{eachROI} ,'sigma ~ ecc  + (1|subject)');


end
lme_temporal
%% compare intercept only vs slope+intercept
lme_temporal_1= fitlme(vertcat(roiTable{:}), 'temporal ~ ecc  + (1|subject)');
lme_temporal_2= fitlme(vertcat(roiTable{:}), 'temporal ~ ecc  + (ecc|subject)');
compare( lme_temporal_1,lme_temporal_2)
% 
for eachROI = 1:length(roiList)
    compare( lme_temporal{eachROI}, lme_temporal2{eachROI})
    
end

%%
fn.fig1 = figure('Renderer', 'painters', 'Position', [100 100 1200 400]);
saveName.fig1 = ['ecc_vs_tw'];


gcolor = turbo(10);
for er = 1:length(roiList)

    subplot(1,length(roiList),er)
    gscatter(roiTable{er}.ecc,roiTable{er}.temporal,roiTable{er}.subject,gcolor,'.',15); hold on
    
    lme_temporal_1= fitlme((roiTable{er}), 'temporal ~ ecc');
    tblnew = table();
    tblnew.ecc = linspace(0,12,12)';
    tblnew.temporal = linspace(4,100,12)';
    [ypred,yCI,DF] = predict(lme_temporal_1,tblnew);
   
    plot(tblnew.ecc,ypred,'color','r','linewidth',3);  hold on

    hold on;
xticks([0 3 6 9 12])
yticks([0 5 10 15 20 25 30])

xticklabels(xticks)
yticklabels(yticks*10);

ylim([0 23]); xlim([0 12]);
if er ==1
ylabel('time window (ms)');
elseif er == 2
    xlabel('eccentrivity (deg)'); 
end
tch_set_axes;
 axis square;
legend off;


    
    
slp  = lme_temporal{er}.fixedEffects;
pval = lme_temporal{er}.Coefficients.pValue;
xloc = 7; yloc = 20;
if pval(2) < 0.05 && pval(2) > 0.01
    text(xloc ,yloc , sprintf('b = %0.2f \np < .05',slp(2)),'FontSize',14);
elseif pval(2) < 0.01 && pval(2) > 0.001
    text(xloc ,yloc , sprintf('b = %0.2f \np < .01',slp(2)),'FontSize',14);
elseif pval(2) < 0.001
    text(xloc ,yloc , sprintf('b = %0.2f \np < .001',slp(2)),'FontSize',14);
else
    text(xloc ,yloc , sprintf('b = %0.2f \np = n.s',slp(2)),'FontSize',14);
end
title(roiList(er))
end



 
end
