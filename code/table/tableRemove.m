function newDT = tableRemove(DT,params)

    for jj=1:width(DT)
        
        paramName   = DT.Properties.VariableNames{jj};
        
        if isfield(params,paramName) %%|| isempty(params.(paramName))
            removerows{jj} = (DT.(paramName) == params.(paramName) );
%             DT(rows,:)=[];
        end
    end
    removerows= removerows(~cellfun('isempty',removerows));
    idx = sum(cell2mat(removerows),2)== size(removerows,2);
    DT(idx,:) = [];
    newDT = DT;

end