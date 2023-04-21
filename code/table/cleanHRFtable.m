function hrfTable = cleanHRFtable(hrfTable)


hrfTable.name = string(hrfTable.name);

hrfTable = hrfTable(~contains(hrfTable.name,'V3d'),:);
hrfTable = hrfTable(~contains(hrfTable.name,'V3v'),:);
hrfTable = hrfTable(~contains(hrfTable.name,'V2d'),:);
hrfTable = hrfTable(~contains(hrfTable.name,'V2v'),:);
hrfTable = hrfTable(~contains(hrfTable.name,'Copy'),:);

% hrfTable = hrfTable(~contains(hrfTable.name,'VO1'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'VO2'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'LO1'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'LO2'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'TO1'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'TO2'),:);
% 
% hrfTable = hrfTable(~contains(hrfTable.name,'IPS0'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'IPS1'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'IPS2'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'IPS3'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'IPS4'),:);
% hrfTable = hrfTable(~contains(hrfTable.name,'IPS5'),:);

mergeTargets = {'IPS','TO','LO','VO'};
for mt = 1:length(mergeTargets)
    hrfTable(contains(hrfTable.name,mergeTargets{mt}),:).name = ...
    repmat(mergeTargets{mt},size(hrfTable(contains(hrfTable.name,mergeTargets{mt}),:).name));
end

hrfTable.name = erase(hrfTable.name,'rh_');
hrfTable.name = erase(hrfTable.name,'lh_');


end