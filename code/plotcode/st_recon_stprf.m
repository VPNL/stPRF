function st_recon_stprf(DT,params,roiList,channel,plotDir)


saveName = sprintf('recon.mat');
saveName = [channel '_' saveName];

if ~isfile(fullfile((plotDir),saveName))
    params.roiName2 = 'all';
    df = tableRemove(DT,params);
    df = defaultThresh(df,0);
    
    params=[]; recon =[];
    id = unique(df.subjID);
    for sj = 1:length(id)
        disp(sj)
        params.subjID =  id(sj);
        [~, ~,stRF] = st_plot_stparam_visualized(df, params,roiList,channel,0,0);
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

end
%%