
function RF = createRF(rf)



x = linspace(0,24,61);
% x = (0.1980:0.1980:20);
x = x - mean(x);
y = x;

% Calculate the spatial sampling parameters
[X,Y] = meshgrid(x,y); XY = [{X},{-Y}];
% Y = -Y;
rf.values = pmGaussian2d(XY{1}, XY{2}, ...
    rf.sigmaMajor,rf.sigmaMinor,rf.Theta, ...
    rf.Centerx0,rf.Centery0);
% rf.values = normMax(rf.values(:));
% RF = reshape(rf.values(:),size(rf.values));
RF = rf.values;
end
