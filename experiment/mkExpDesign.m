function ExpDesign = mkExpDesign(nSpace,nTemporal)
% counterBalance spatial and temporal conditions.
% Tries to minimize presentation order bias
% nSpace = 36; nTemporal = 9;

% output 
% ExpDesign: temporal sequence of [Run X Spatial Loc]
%
%
[space, time] = BalanceFactors(1, 1, 1:nSpace,1:nTemporal);

for steps = 1:36
    ord(:,steps) = time(find(space==steps));
end

minimize = std(sum(ord,2));
while minimize > 4
    [space, time] = BalanceFactors(1, 1, 1:36,1:9);

for steps = 1:36
    ord(:,steps) = time(find(space==steps));
end

minimize=std(sum(ord,2));
ExpDesign = ord;

end

