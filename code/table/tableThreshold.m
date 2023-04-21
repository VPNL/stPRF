
function newDT = tableThreshold(DT,cutoff,level)

if notDefined('level')
    level = 'lower' ; % remove everything lower this
end

for jj=1:width(DT)
    
    paramName   = DT.Properties.VariableNames{jj};
    
    if isfield(cutoff,paramName) %%|| isempty(params.(paramName))
        if strcmp(level,'lower')
            
            if iscell(DT.(paramName))
                rows = (cell2mat(DT.(paramName)) < cutoff.(paramName) );
                rows = sum(rows,2) < 1;
            else
                rows = (DT.(paramName) < cutoff.(paramName) );
                rows = sum(rows,2) < 1;
            end
            
            fprintf('[%s]: removing voxels below %d %s \n',mfilename, cutoff.(paramName),paramName)
            
            
        elseif strcmp(level,'upper')
            if iscell(DT.(paramName))
                rows = (cell2mat(DT.(paramName)) > cutoff.(paramName) );
                rows = sum(rows,2) < 1;
            else
                rows = (DT.(paramName) > cutoff.(paramName) );
                rows = sum(rows,2) < 1;
                
            end
            %             rows = (DT.(paramName) > cutoff.(paramName) );
            fprintf('[%s]: removing voxels above %d %s \n',mfilename, cutoff.(paramName),paramName)
            
        end
        
        DT = DT(rows,:);
        fprintf('[%s]: %s thresholded %d voxels\n',mfilename,paramName,sum(~rows))
        
    end
    newDT = DT;
end

end