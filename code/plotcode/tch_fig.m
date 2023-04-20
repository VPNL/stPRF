function fig = tch_fig(name, pos)

if nargin < 1; name = ''; end
if nargin < 2; pos = [.1 .1 .5 .5]; end
% pos = rectify(pos);
fig = figure('Name', name, 'Color', 'w', 'Units', 'centimeters', 'Position', pos);
tch_set_axes; hold on;

end