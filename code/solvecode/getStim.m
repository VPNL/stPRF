function  stim = getStim(params,doDown)

if notDefined('doDown')
    doDown = 1;
end


%% downsample stimulus and reshape according to the temporal sampling rate
if doDown == 1
    downSampleRate = 1000/params.analysis.temporal.param.fs;
    stim = params.stim(1).images_unconvolved';
    stim_all = zeros(size(stim,2),size(stim,1)/10,length(params.stim));
    if params.analysis.temporal.param.fs ~= 1000
        for es = 1:length(params.stim)
%             stim = full(logical( params.stim(es).images_unconvolved'));
            stim = full(( params.stim(es).images_unconvolved'));
            stim_all(:,:,es) = single(downsample(stim,downSampleRate)'); %3721 21000 6
        end
    end

else
    stim_all = [];
    for es = 1:length(params.stim)
%         stim_all(:,:,es) = full(logical( params.stim(es).images_unconvolved));
        stim_all(:,:,es) = full(( params.stim(es).images_unconvolved));

    end

end

%% check gpu
if params.useGPU
    g = gpuDevice();
    reset(g);
    if g.DeviceSelected %canUseGPU()
        stim = gpuArray(stim_all);
    end
else
    stim = stim_all;
end

end