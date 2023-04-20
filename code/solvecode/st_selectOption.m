% Selects one of default option bundles
% current working version is analysis option 2 or 97

function optionBundle = st_selectOption(analysisoption)


switch analysisoption
    
    % Maximal Version
    % no spaital smoothing
    % currentWorking Version
    case 1
        rmName = 'GLMDenoised';
        numberStimulusGridPoints = 30;
        keepAllPoints = true;
        coarsetofine = true;     % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 1;       
        
    % spaital smoothing
    case {2, 97}
        rmName = 'MotionComp_RefScan1';
        numberStimulusGridPoints = 30;
        keepAllPoints = false;
        coarsetofine = true;     % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 1;              % do cross validation;
        normrf = 0;
  
        % without crossvalidation
    case 3
        rmName = 'MotionComp_RefScan1';
        numberStimulusGridPoints = 30;
        keepAllPoints = false;
        coarsetofine = false;    % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 0;              % do cross validation;
        normrf = 0;

    case 4
        rmName = 'MotionComp_RefScan1';
        numberStimulusGridPoints = 30;
        keepAllPoints = true;
        coarsetofine = false;    % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 1;              % do cross validation;
        normrf = 1;
        
    case 5 
        rmName = 'MotionComp_RefScan1';
        numberStimulusGridPoints = 30;
        keepAllPoints = true;
        coarsetofine = false;     % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 0;              % do cross validation;

    case 6  % GLMdenoise and smooth
        rmName = 'GLMDenoised';
        numberStimulusGridPoints = 30;
        keepAllPoints = true;
        coarsetofine = true;     % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 0;              % do cross validation;
        
    case 7  % GLMdenoise and no smoothing
        rmName = 'GLMDenoised';
        numberStimulusGridPoints = 30;
        keepAllPoints = true;
        coarsetofine = false;     % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 0;              % do cross validation;

    case 8  %   % 3-fold cv
        rmName = 'GLMDenoised';
        numberStimulusGridPoints = 30;
        keepAllPoints = true;
        coarsetofine = false;     % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 1;              % do cross validation;

    case 9  %   % HRF version3-use optimized HRF
        rmName = 'GLMDenoised';
        numberStimulusGridPoints = 30;
        keepAllPoints = true;
        coarsetofine = false;     % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        hrfver = 3;
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 0;              % do cross validation;

    case 10  %   % HRF version3-use optimized HRF
        rmName = 'GLMDenoised';
        numberStimulusGridPoints = 30;
        keepAllPoints = true;
        coarsetofine = false;     % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        hrfver = 3;
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 1;              % do cross validation;

    % in-development
    % allAvged TC analysis
    case 98 % development
        rmName = 'allAvg';
        numberStimulusGridPoints = 20;
        keepAllPoints = false;
        coarsetofine = false;    % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 0;              % do cross validation;
        
        
    % Quick check if the code works
    % 1) whether grid can be created
    % 2) whether cross validation works
    case 99
        
        rmName = 'MotionComp_RefScan1';
        numberStimulusGridPoints = 1;
        keepAllPoints = false;
        coarsetofine = false;    % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 1;              % do cross validation;
        
        
    % Quick check if the code works
    % 1) whether grid can be created
    % 2) whether cross validation works
    % 3) Spatial Blur
    case 100
        rmName = 'MotionComp_RefScan1';
        numberStimulusGridPoints = 1;
        keepAllPoints = false;
        coarsetofine = true;    % spatial smoothing
        coarsetblurparams = 0;   % temporal smoothing (amount to decimate)
        calcPC = 1;              % covert to percent-signal change
        cache  = 1;              % save grid outputs
        cv     = 1;              % do cross validation;
        
    otherwise
        msg = sprintf('[%s] select valid option',mfilename);
        error(msg);
end


optionBundle.rmName                   = rmName;
optionBundle.numberStimulusGridPoints = numberStimulusGridPoints;
optionBundle.keepAllPoints            = keepAllPoints;
optionBundle.coarsetofine             = coarsetofine;
optionBundle.coarsetblurparams        = coarsetblurparams;
optionBundle.calcPC                   = calcPC;
optionBundle.cache                    = cache;
optionBundle.cv                       = cv;

if notDefined('hrfver')
    hrfver = [];
end
optionBundle.hrfver                   = hrfver;

% optionBundle.normrf                   = normrf;

end
