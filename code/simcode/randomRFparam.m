function RFparam = randomRFparam(StimSize,nVoxels)
rng('default') % For reproducibility

RFxy = normrnd(0,StimSize/4,[nVoxels,2]);

% if we use linlog (default) we can specify where the cutoff lies. This
RFsigma = boundUniformRand([0.2 StimSize/4],nVoxels);

RFparam = [ round(RFxy,2)  round(RFsigma,2)];
end
