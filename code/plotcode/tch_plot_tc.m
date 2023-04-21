function [h, ymin, ymax] = tch_plot_tc(x, y, lw, lcol, ecol, err_type)
% Generates clean version of a simple line plot with error region
% 
% INPUTS
%   1) x: vector of x-axis values
%   2) y: matrix of data (N x P)
%   3) lw: line width (default = 1)
%   4) lcol: plotting color for line
%   5) ecol: plotting color for error region
%   6) err_type: type of error region ('std' or 'sem')
% 
% OUTPUTS
%   1) h: handle to line object
%   2) ymin: minimum value of bar - err
%   3) ymax: maximum value of bar + err
% 
% AS 5/2017

% check inputs
if nargin < 3 || isempty(lw); lw = 2; end
if nargin < 4 || isempty(lcol); lcol = [0 0 0]; end
if nargin < 5 || isempty(ecol); ecol = [.5 .5 .5]; end
if nargin < 6 || isempty(err_type) || size(y, 1) < 2
    err_flag = 0; err_type = 'std';
else
    err_flag = 1;
end

% compute mean and error of y
y_mean = nanmean(y, 1);
if nargin > 4 && size(y, 1) > 1
    err = std(y);
    if strcmpi(err_type, 'sem')
        err = err / (sqrt(size(y, 1) - 1));
    end
else
    err = zeros(1, size(y, 2));
end

% plot line and error region
if err_flag
    xvec = [x fliplr(x)]; yvec = [y_mean + err fliplr(y_mean - err)];
    me = patch(xvec, yvec, ecol, 'LineStyle', 'none'); alpha(me, 0.2);
    me.Annotation.LegendInformation.IconDisplayStyle = 'off';

end
hold on; ma = plot(x, y_mean,'Color', lcol, 'LineWidth', lw);
ymin = min(y_mean - err); ymax = max(y_mean + err);
if err_flag == 1;
    h = ma;
else
    h = ma;
end

end
