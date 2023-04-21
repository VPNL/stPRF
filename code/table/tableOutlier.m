
function newDT = tableOutlier(DT,targetParam,outMethod)

if notDefined('outMethod')
   outMethod =  'median';
end
for jj=1:width(DT)
    
    paramName   = DT.Properties.VariableNames{jj};
    
    if strcmp(targetParam,paramName)%%|| isempty(params.(paramName))
        if iscell(DT.(paramName))
            rows = ~isoutlier(cell2mat(DT.(paramName)),outMethod);
        else
            rows = ~isoutlier(DT.(paramName),outMethod);
        end
        rows = (sum(rows,2) == size(rows,2));
        DT = DT(rows,:);
        fprintf('[%s]: removing %d %s outliers\n',mfilename,sum(~rows),paramName)
    end
    newDT = DT;
end

end