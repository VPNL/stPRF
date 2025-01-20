function [trial,nStim] = st_mkTrial(stimprm,fs)

% function [stimulus, time] = st_fMRIStimulus(with0OrNot)

% INPUT -------------------------------------------
% stimprm: [variations X prms]
%     prms = [on off trialDur framerate]
%     on,        unit: nFrame
%     off (isi), unit: nFrame
%     trialDur,  unit: second
%     framerate, unit: HZ
%     


% OUTPUT -------------------------------------------

% trial : nTrial X time (ms)

% History :

%% Example


%% Assess inputs

% possibleInputs = {'with0', 'without0'};
% assert(ismember(lower(with0OrNot), possibleInputs), 'Un-identifiable input.');

%% Pre-defined variables

% fs = 1000; %1KHz sample rate
nStimulus = size(stimprm, 1);


%% compute stimulus
trial = cell(nStimulus,1);

for istim = 1 : nStimulus
    prm = stimprm(istim,:);
    
    trial_dur = prm(3);
    stim_on = prm(1)/ prm(4);
    stim_off = prm(2)/ prm(4);
    stim_num =  trial_dur / ( stim_on +  stim_off);
    hz =  stim_num/ trial_dur;
    
    
    t = 0:1/fs:trial_dur; % time vector
    w =  stim_on;
    d = w/2:1/ hz:trial_dur ;% Delay starts at t = 0 and
    each_stim = pulstran(t,d,'rectpuls',w);
    each_stim = each_stim(1:end-1); % cut off the last bit
    
    if  hz == inf
        each_stim = ones(size(each_stim));
    end
    
    trial{istim} = each_stim;
    nStim{istim} = stim_num;    
end


end