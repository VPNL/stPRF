function [fn,saveName] = st_plot_param_scatter(DT, params,roiList,targetVar)

rng('default')

sampleSize = 50;

cc = getmycolors;

% initial filter
df = tableFilter(DT,params);

% give number to ROIs
df = st_tablegiveROInumber(df, roiList);

id = unique(df.subjID);  


if strcmp(targetVar,'sigma')
    x1 = linspace(0,8);
elseif strcmp(targetVar,'ecc')
    x1 = linspace(0,12);
elseif strcmp(targetVar,'phase')
    x1 = linspace(-3.14,3.14);
else
    x1 = linspace(-12,12);
end
  
%%
counter =1; 
for er = 1:length(roiList)
    params = [];
    params.roiNumber =  er;
    df_temp = tableFilter(df,params);
    df_temp = convertvars(df_temp,{'sigma','x0','y0','ecc'},'double');


    for sj = 1:length(id)
        params = [];
        params.subjID =  id(sj);
        params.modelNumber = 1;
        df1 = tableFilter(df_temp,params);
        params.modelNumber = 2;
        df2 = tableFilter(df_temp,params);
        params.modelNumber = 3;
        df3 = tableFilter(df_temp,params);

        if size(df1.(targetVar)) < sampleSize
            sampdata1 = [df1.(targetVar),df2.(targetVar)];
            sampdata2 = [df1.(targetVar),df3.(targetVar)];
        else
            sampdata1 = datasample([df1.(targetVar),df2.(targetVar)],sampleSize,'Replace',true);
            sampdata2 = datasample([df1.(targetVar),df3.(targetVar)],sampleSize,'Replace',true);

        end
        
        p1=polyfit(df1.(targetVar),df2.(targetVar),1);
        fitSlopes{sj,er}  = p1(1);
        fitData1{sj,er} = polyval(p1,x1);
        boxData1{sj,er} = sampdata1;

        p2=polyfit(df1.(targetVar),df3.(targetVar),1);
        fitData2{sj,er} = polyval(p2,x1);
        boxData2{sj,er} = sampdata2; %[df1.sigma df2.sigma ];
        counter= counter+1;
    end 
end





%%
fn=figure(); clf;
ScrSz = get(0, 'screensize');
set(gcf, 'position',  [ScrSz(1) ScrSz(2) ScrSz(3)/2 ScrSz(4)/3]);
saveName = ['scatter-' targetVar];

pos=[];
pos{1} = [1 2 7 8];
pos{2} = 3; pos{3} = 9; pos{4} = 4;
pos{5} = 10; pos{6} = 5; pos{7} = 11;
pos{8} = 6; pos{9} = 12; 

for er = 1:length(roiList)
    
 
         ax =   subplot(2,6,pos{er});

bd = cell2mat(boxData1(:,er));
scatter(bd(:,1),bd(:,2),10, ...
    'MarkerFaceColor',cc(2,:),...
    'MarkerFaceAlpha',.5, ...
    'MarkerEdgeColor','none',...
    'MarkerEdgeAlpha',.1); hold on



bd = cell2mat(boxData2(:,er));
scatter(bd(:,1),bd(:,2),10, ...
    'MarkerFaceColor',cc(3,:),...
    'MarkerFaceAlpha',.5, ...
    'MarkerEdgeColor','none',...
    'MarkerEdgeAlpha',.1); hold on


plotData = cell2mat(fitData1(:,er));
plotData( ~any(plotData,2), : ) = [];  %rows
tch_plot_tc(x1,plotData,1.5,cc(2,:),cc(2,:),'sem'); hold on;


plotData = cell2mat(fitData2(:,er));
plotData( ~any(plotData,2), : ) = [];  %rows
plot(linspace(-12,12,10),linspace(-12,12,10),'k--','linewidth',2); hold on;
tch_plot_tc(x1,plotData,1.5,cc(3,:),cc(3,:),'sem'); hold on;

yticks([min(x1) mean(x1) max(x1)]);
xticks([min(x1) mean(x1) max(x1)]);
if strcmp(targetVar,'phase')
    xticklabels({'-\pi','0','\pi'})
    yticklabels({'\pi','0','\pi'})
end
xlim([min(x1) max(x1)]); ylim([min(x1) max(x1)]);

axis square;
tch_set_axes;
title(roiList(er),'Fontsize',14);

if er > 1
    set(ax,'fontsize',14)
end
end


end
