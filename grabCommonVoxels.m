function dff = grabCommonVoxels(DT,targetTmodel)

% need 3 models for now
% grab common voxels across all three models and filter it out

df = []; params=[];
sj = unique(DT.subjID);
for tm = 1:3
    params.tmodel = targetTmodel{tm};
    df{tm} = tableFilter(DT,params); 
end

%%
tmp=[];
for tm = 1:3
    tmp{tm} = table(df{tm}.coords, df{tm}.subjID , df{tm}.roiName);
end
inter = intersect(intersect(tmp{1},tmp{2},'stable'),tmp{3},'stable');

for tm = 1:3
    [~,ai,~]=intersect(tmp{tm},inter,'stable');
    dff{tm} = df{tm}(ai,:);
end


fprintf('[%s]: grabbing common voxels \n',mfilename)
for tm= 1:3
    fprintf('[%s]: removed %d from %s \n',mfilename,(size(df{tm},1) - size(dff{tm},1)),targetTmodel{tm})
end

dff = [dff{1}; dff{2}; dff{3}];



end