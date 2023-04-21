function st_recon_stprf_ecc(DT,params,roiList,eccBand,plotDir)

saveName = sprintf('%s_recon.mat',eccBand);

if ~isfile(fullfile(fileparts(plotDir),'recon',saveName))
    params.roiName2 = 'all';
    df = tableRemove(DT,params);
    df = defaultThresh(df,0);
    
    params=[]; recon =[];
    id = unique(df.subjID);
    for sj = 1:length(id)
        disp(sj)
        params.subjID =  id(sj);
        [~, ~,stRF] = st_plot_stparam_visualized(df, params,roiList,'sustained',0,0);
        recon{sj} = stRF;
    end
    
    
    id = unique(DT.subjID);
    
    aRF=[];
    for sj = 1:length(id)
        tmp = cellfun(@(x) mean(x,3),recon{sj},'UniformOutput',false)';
        aRF = [aRF tmp];
    end
    
    avgrecon = [];
    for er = 1:length(roiList)
        tmp = aRF(er,:);
        Y = cat(3,tmp{:});
        data = nanmean(Y,3);
        
        % marginalize
        data = (normMax(sum(data,2)) * normMax(sum(data))) ;
        
        avgrecon{er} = data;
        
    end
       
    save(fullfile(plotDir,saveName),'recon','avgrecon','-v7.3')
end

%% plot it
% close all; clc;
saveNames =   {  'fov_recon.mat', 'peri_recon.mat'};

counter = 1;
fn = figure('position',[         0         0        900     500])

for ee = 1:length(saveNames)
    saveName = [saveNames{ee}];
    load(fullfile((plotDir),saveName),'avgrecon')
    
    id = unique(DT.subjID);
    
   
    
    for eachROI = 1 :length(out)
        ax{eachROI} = subplot(2,3,counter);


        
        
        data = avgrecon{eachROI};
        imagesc(1:501,y,data); hold on;
        contour(1:501,y,data,[0.1,0.5,0.9],'k'); hold on;
        set(gca, 'YDir', 'normal'); 
        
        axis square; axis tight;
        tch_set_axes;
        set(gca,'linewidth',2)
        
        if counter < 4
            title(roiList{eachROI})
        end
        box off;
        
        colormap(ax{eachROI},bluewhitered);
        
%         cbh = colorbar('h');
%         cbh.Limits;


            grid on; %grid minor;
            aa = gca;
            set(gca,'GridLineStyle','-')
%             set(gca,'MinorGridLineStyle','-')

%             aa.YRuler.MajorTickValues = [0:20:500]; %just like major ticks

%             aa.YRuler.MinorTickValues = [0:20:500]; %just like major ticks
%             aa.XRuler.MinorTickValues = [0:20:500]; %just like major ticks

%             aa.YRuler.MinorTickValuesMode = 'auto'; %or 'manual'

%             set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color

%             aa.GridColor = [1, 1, 1];  % [R, G, B]

%         xticks(0:10:500);  xticklabels({'0','100','200','300', '400'});
       
        xticks(0:10:500); xticklabels({'0','100','200','300', '400'});
        set(gca,'FontName','Times','fontsize',17)

        xlim([0 35])    
        ylim([-5 5])
        counter = counter+1;
    end


    
end
% [~,n,~]=fileparts(saveName);
% printnice(fn,[1 300],plotDir,n)

end
%%
