function oct = draw_oct(x, y, w, h, fill_col, line_col)
% Plots an octogon centered at (x, y) of width 2w and height 2h
% INPUTS
%   1) x: x-axis value of octogon origin
%   2) y: y-axis value of octogon origin
%   3) w: width from origin to horizontal extremes
%   4) h: height from origin to vertical extremes
%   5) fill_col: octogon fill color
%   6) line_col: outline color
% 
% OUTPUT
%   1) oct: handle to octogon fill object
% 
% AS 6/2017

t = (1/16:1/8:1)' * 2 * pi;
xx = w * sin(t) + x;
yy = h * cos(t) + y;
oct = fill(xx, yy, fill_col, 'EdgeColor', line_col);

end
