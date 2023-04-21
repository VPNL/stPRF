function [fn, savename] = st_plot_optHRF(hrfTable,sumTable)

% hrfTable = cleanHRFtable(hrfTable);
% hrfTable= tableRemoveHRFDuplicateCoords(hrfTable);
% sumTable= tableSumHRF(hrfTable);


%%
fn = figure('position',[0 0 1500 500]);
savename = ['hrfs'];

% figure panel 1
params=[];
roiList = {'V1','V2','V3','hV4','VO','LO','TO','V3AB','IPS'};
gc = (getmycolors(3));
fs =10;
subplot(2,5,[1 2 6 7])
for r = 1:1:length(roiList)
    roi = roiList(r);
    params.name=roi;
    newDT = tableFilter(sumTable,params);
    hrfs = normMax(newDT.hrfs')';
    tch_plot_tc(1:201,hrfs,2,gc(r,:),gc(r,:),'std'); hold on;
end

vistaParams     = [5.4 5.2 10.8 7.35 0.35];
tSteps = 0:1/fs:20;
values = rmHrfTwogammas(tSteps, vistaParams);
plot(values,'k--','linewidth',2);

newLegend = roiList;
newLegend{end+1} = 'Vistasoft HRF';
legend(newLegend); legend box off;
axis square;

hrfaxes();


% figure panel 2
roiList = {'V1','LO'};
gc = jet(10);
mapping = [3,4];
for r = 1:length(roiList)
    
    
    subplot(2,5,mapping(r))

    for s = 1:10
        
        roi = roiList(r);
        params.subj=s;
        params.name=roi;
        newDT = tableFilter(hrfTable,params);
        tch_plot_tc_norm(1:201,newDT.hrfs,2,gc(s,:),[.7 .7 .7],'std'); hold on;
        hrfaxes();
        
    end
    
    text(130,0.8,roi,'fontweight','bold', 'FontSize', 14 ...
        , 'FontName', 'Helvetica');
    tch_set_axes

end


% panel 3
params=[];
gc = (getmycolors(3));
rng('default');
roiList = {'V1','V2','V3','hV4','VO','LO','TO','V3AB','IPS'};

for s = 3
    for r =[1, 6]
        if r == 1
            subplot(2,5,8)
        elseif r ==6
            subplot(2,5,9)
        end
        roi = roiList(r);
        params.subj=s;
        params.name=roi;
        newDT = tableFilter(hrfTable,params);
        
        for t=  randperm(height(newDT),500)
            plot(newDT.hrfs(t,:),'color',[gc(r,:) 0.1],'linewidth',0.2) ; hold on;
        end
        text(130,0.8,roi,'fontweight','bold', 'FontSize', 14 ...
            , 'FontName', 'Helvetica');
        tch_plot_tc_norm(1:201,newDT.hrfs,3.5,[1 1 1]); hold on;
        tch_plot_tc_norm(1:201,newDT.hrfs,2,[0 0 0]); hold on;
        hrfaxes();
    end
    
    
end

end
%%