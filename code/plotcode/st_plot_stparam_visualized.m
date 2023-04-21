function [fn,saveName,all_stRF] = st_plot_stparam_visualized(DT, params, roiList,channel, normVersion,applyExp)
% params=[]; close all;
fn = figure('Renderer', 'painters', 'Position', [0 0 1800 500]);
tms1 = {'\tau_1', '\itw', '\tau_2', '\itn_D_N', '\sigma_D_N'};
tms2 = {'\tau', '\itn'};

% plot info
polar_radius = 12;
polar_dim = 61;
x = linspace(0,polar_radius*2,polar_dim);
x = x - mean(x); y = x;
tWindow = 501;

% initial filter
DT = tableFilter(DT,params);

% give number to ROIs
DT = st_tablegiveROInumber(DT, roiList);

all_stRF = cell(size(unique(DT.subjID),1),size(roiList,2));
for targetModel = 3
    params.modelNumber = targetModel;
    df = tableFilter(DT,params);
    
    if targetModel == 3
        thresh = [];
        thresh.temporal = 99;
        df = tableThreshold(df,thresh,'upper'); % remove anything upper;
        
        thresh = [];
        thresh.temporal = 4;  % max ecc
        df = tableThreshold(df,thresh,'lower'); 
    end
    
    % % unpack table cells to arrays
    df = tableArrayFormat(df);
    
    for eachROI = 1:length(roiList)
        subplot_tight(1,length(roiList),eachROI,[0.1,0.03])
        params.roiNumber = eachROI;
        df2 = tableFilter(df,params);
        if ~isempty(df2)
            df2 = st_reconstruct(df2,1,0);
        end
        %%
        
        stRF = zeros( polar_dim, tWindow, height(df2), 'single' );
        
        for erf = 1:height(df2)
            cRF = cell2mat(df2(erf,:).sRF);
            IRF = cell2mat(df2(erf,:).tRF);
            
            if strcmp(channel,'sustained')
                IRF = IRF(:,1)';
            elseif strcmp(channel,'transient')
                IRF = IRF(:,2)';
            end
            
            [ ~, ix ] = max(cRF(:));
            [ i1,i2] = ind2sub( size(cRF), ix );
            px = cRF(i1,:)';
            py = cRF(:,i2);
            
            if normVersion ==1
                st = (py)*(IRF);
                st = normMax(st(:));
                stRF(:,:,erf) = reshape(st,size(py,1),tWindow);
            else
                stRF(:,:,erf) =((py)*(IRF));
            end
            
            
            if applyExp == 1
                st = ((py)*(IRF));
                st = real(power( st,df2(erf,:).exponent));
                stRF(:,:,erf) = st;
            end
            
            
        end
        p_stRF= mean(stRF,3);
        
        
        p_stRF = p_stRF./max(max(p_stRF));
        all_stRF{eachROI} = stRF;
        
        
        imagesc(1:size(sum(stRF,3),2),y,p_stRF); colormap(bluewhitered);  hold on;
        % contour(1:size(sum(stRF,3),2),y,p_stRF,[0.05 0.25 0.5 0.75],'k'); hold on;
        contour(1:size(sum(stRF,3),2),y,p_stRF,5,'k'); hold on;
        
        set(gca, 'YDir', 'normal'); grid on;
        axis square; axis tight;
        xticks(0:20:500);  xticklabels({'0','200','400','600'});
        
        tch_set_axes;
        xlim([0 60])
        ylim([-12 12])
        set(gca,'linewidth',2)
        title(roiList{eachROI})
    end
    
    
    
end

tch_set_axes;
legend(roiList);
legend  off;
saveName = ['stparams-vis.png'];


end

