
function newDT = tableFilter(DT,params)
    for jj=1:width(DT)

        paramName   = DT.Properties.VariableNames{jj};

        if isfield(params,paramName) %%|| isempty(params.(paramName))
            rows = (DT.(paramName) == params.(paramName) );
            DT = DT(rows,:);
        end
        newDT = DT;
    end
end