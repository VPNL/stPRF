function params = displayParams
% For Scanner NEC LCD2080UX display
% Driven directly by Powerbook laptop via a labelled DVI-VGA cable, a VGA F-F adapter,
% a Y-shape VGA splitter, with the other end hang empty, a VGA-5 BNC converter, 5 BNC
% connectors at control room<->scanner room panel, 50ft BNC cables, and 5 BNC-VGA converter.
% LCD params: Brightness = 100, Contrast = 50, colorTemp = NATIVE
% 800x600 60Hz, for both Apple's panel and LCD. Others in default.
% Last calibrated: 09-11-2003. Using USB-serial UC232A, PR650 remote at 2deg mode (center black dot)
% PR650 about 280 cm away, center at screen.

% Critical parameters
params.numPixels = [1920 1080];
params.dimensions = [39 22];
params.distance = 60;%backbore
params.frameRate = 60;
params.cmapDepth = 10;
params.screenNumber = 1;

% Descriptive parameters
params.monitor = 'propixx';
params.position = '3T scanner';





% parameters which can make programming the display easier
% params.flipLR = 1;
% params.flipUD = 0;

return

