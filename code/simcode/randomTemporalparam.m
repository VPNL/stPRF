function tParam = randomTemporalparam(temporalModel,nVoxels)
rng('default') % For reproducibility

switch temporalModel
    case 'spatial'
        searchRange =[];
    case 'CST'
        % 2 extra params to solve:
        % [temporal delay 2) exponent]
        searchRange = [
            4   0.1; ...
            100  1];
    case 'DN-ST'
        %  ["tau1", weight, "tau2", "n", "delay/sigma"]
        searchRange = [
            0.01  0  0.01   1  0.01; ...
            1     0     1   6   0.5];
end


tRange = [min(searchRange); max(searchRange)];
nParam = size(tRange,2);
for ep = 1:nParam
    tParam(:,ep)= round(boundUniformRand(tRange(:,ep),nVoxels),2);
end


switch temporalModel
    case 'CST'
        tParam = [tParam(:,1) tParam]; % we want same time-to-peak (tau) for both sustaiend and transient
end



end

