function [pred,beta,R2]=stReFit(params,stim,data,t)
%%
    fracAlpha = params.analysis.optim.ridgeAlpha;
    nChan = getChanNumber(params);

    predictions = stPredictBOLDFromStim(params, stim);
    predBOLD = predictions.predBOLD;
    
    if size(stim,3) ==1
        predBOLD = squeeze((predBOLD));
    else
        predBOLD = squeeze(concatRuns(predBOLD));
    end
    predBOLD = [predBOLD t.trends];
    [betas,~,~] = fracridge(gather(predBOLD),fracAlpha,gather(data),[],1);
    
    pred = gather(predBOLD * betas);
    beta = gather(betas(1:nChan));
    R2 = calccod(gather(pred),gather(data))/100;
end

