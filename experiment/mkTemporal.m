function [img_sequence, nimg, img_durations]= mkTemporal(ExpDesign,temporal)
% counterBalance spatial and temporal conditions.
% Tries to minimize presentation order bias
% nSpace = 36; nTemporal = 9;

% output 
% ExpDesign: temporal sequence of [Run X Spatial Loc]
%
%

%%


% nFrames
ifi = 1/temporal.fs;

if ifi ==60
    ons = [8 8  8  2 12 52 16 48 300]';
    offs= [2 12 52 8 8  8  4  12   0]';
elseif ifi == 75
    ons =  [10, 10, 10, 5 , 15, 65,    17, 50, 375]';
    offs = [5,  15, 65, 10, 10, 10,     8, 25, 0]';
end
trial_dur = temporal.trial_dur;
stimprm = zeros(length(ons),4);
for i = 1:length(ons)
    stimprm(i,1)=ons(i);
    stimprm(i,2)=offs(i);
    stimprm(i,3)=trial_dur;
    stimprm(i,4)=ifi;
end


[trial,nstim]=st_mkTrial(stimprm,1/temporal.fs);
% [trial,nstim]=st_mkTrial(stimprm,1000);


for tp = 1:size(stimprm,1)
    temp{tp} = [ones(1,stimprm(tp,1)) zeros(1,stimprm(tp,2))];
end

img_sequence=[];
img_durations =[];
for ii = 1:length(ExpDesign)
    img_sequence{ii} = temp{ExpDesign(ii)};
    img_durations{ii} = repmat(round(stimprm(ExpDesign(ii),1:2)*temporal.fs*1000)',1,single(nstim{ExpDesign(ii)}));
    nimg{ii}        = single(nstim{ExpDesign(ii)});
end

img_durations = cell2mat(img_durations)';

end

