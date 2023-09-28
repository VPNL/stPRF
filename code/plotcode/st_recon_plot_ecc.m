function [fn,savename]=st_recon_plot_ecc(DT,params,roiList,plotDir)



polar_radius = 12;
polar_dim = 61;
x = linspace(0,polar_radius*2,polar_dim);
x = x - mean(x); y = x;

% close all; clc;
saveNames =   {  'fov_recon.mat', 'peri_recon.mat'};

counter = 1;
fn = figure('position',[         0         0        900     500]);
savename = 'fov_vs_peri';


for ee = 1:length(saveNames)
    saveName = [saveNames{ee}];
    load(fullfile((plotDir),saveName),'avgrecon')
    
    id = unique(DT.subjID);
    
   
    
    for eachROI = 1:3 %length(avgrecon)
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
        
        grid on; %grid minor;
        aa = gca;
        set(gca,'GridLineStyle','-')

        xticks(0:10:500); xticklabels({'0','100','200','300', '400'});
        set(gca,'FontName','Times','fontsize',17)

        xlim([0 35])    
        ylim([-5 5])
        counter = counter+1;
    end
    
   
end

end
%%
