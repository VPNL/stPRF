function params = displayParams
% Critical parameters
% for projector
params.numPixels = [1280 1024];
params.dimensions = [38 20]; % [46 20] 
% params.dimensions = [46 20] ; % [46 20] 

params.distance = 51; % 32 channel -> 43 cm // 16 channel -> 25 cm
params.frameRate = 60;
params.cmapDepth = 10;
params.screenNumber = 1;

% Descriptive parameters
params.monitor = 'CNI_projector';
params.position = '3T scanner';
params.flipLR = 1;


% % for LCD
% params.numPixels = [1920 1080];
% params.dimensions = [104 59];
% params.distance = 277; % 32 channel -> 43 cm // 16 channel -> 25 cm
% params.frameRate = 60;
% params.cmapDepth = 10;
% params.screenNumber = 1;

% % Descriptive parameters
% params.monitor = 'CNI_LCD';
% params.position = '3T scanner';
% 




% parameters which can make programming the display easier
% params.flipUD = 0;

return

