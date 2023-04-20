%  unpack table cells to arrays
%  KIS 2021 
%%
function newT = tableArrayFormat(df)
targetVariables = df.Properties.VariableNames;

targetVariables = setdiff(targetVariables, {'sRF', 'tRF'});

for vn = df.Properties.VariableNames
    paramName=vn;
    if max(ismember(targetVariables,paramName))
        if iscell(df.(paramName{:}))
            newT.(paramName{:}) = cell2mat(df.(paramName{:}));
        else
            newT.(paramName{:}) = df.(paramName{:});
        end
    end
    
end
newT=struct2table(newT);

end