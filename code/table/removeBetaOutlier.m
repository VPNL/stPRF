function df = removeBetaOutlier(DT)

targetTmodel =  {'1ch_glm','1ch_dcts', '3ch_stLN'};
for tm = 1:3
    params.tmodel = targetTmodel{tm};
    df{tm} = tableFilter(DT,params);
    df{tm} = tableOutlier(df{tm},'beta');
end
df = [df{1}; df{2}; df{3}];

end