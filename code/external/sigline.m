function sigline(xs,lbl,h,yv)
%SIGLINE plots a line of statistical significance on the current axis
%   sigline(xs) plots a significance line between the value in the 2D
%   vector xs along the x-coordinate. 
% 
%   sigline(xs,lbl) places a text lbl beside the significance line.
% 
%   sigline(xs,lbl,h) attempts to plot the significance line above a
%   highest value in a graphics object with handle h. If h is not a
%   graphics handle, the line will be plotted above the value of h. 
% 
%   sigline(xs,...,vy) give a way to specify that the line be plotted above
%   vy which should not be assumed to be a graphics object. This is
%   optional and overrides h; but probaly won't ever be needed. The vy is a
%   value on the y-axis above which the sig line is plotted. It is only
%   useful for a possibly, rare occurance of a given y in h unwantedly
%   matching a graphic object in value.
% 
%   Examples #1: Plot 5 random points and add a significance line between
%   point 2 and 4.
%   plot(rand(5,1))
%   sigline([2,4])
% 
%   Examples #3: Plot 5 random points as bar chart and add a significance
%   line between point 2 and 3 and print p=0.05 beside it.
%   bar(rand(5,1))
%   sigline([2,3],'p=0.05')
% 
%   Examples #4: Plot 5 random points and add a significance line between
%   point 2 and 4 and specify the handle to be used to obtaine the value
%   above which the line will be plotted.
%   h=plot(rand(5,1))
%   sigline([2,4],[],h)
% 
%   Examples #5: Plot 5 random points and add a significance line between
%   point 2 and 4 and a value above which the line will be plotted.
%   plot(rand(5,1))
%   sigline([2,4],[],0.9)
%
%   TODO: To avoid messing up legend, we should redo this function using
%   drawing rather than plots? A work aroound for now is to plot legends b4
%   calling this function.
if nargin==1
    y=gety(gca);
    lbl=[];
elseif nargin==2
    y=gety(gca);
elseif nargin==3
    y=gety(h);
elseif nargin==4
    y=yv;
end
% Now plot the sig line on the current axis
hold on
xs=[xs(1);xs(2)];
plot(xs,[1;1]*y*1.1,'-k', 'LineWidth',1);%line
if lbl
    text(mean(xs), y*1.15, lbl,'HorizontalAlignment', 'center')% the sig star sign
else
    plot(mean(xs), y*1.15, '*k')% the sig star sign
end 
plot([1;1]*xs(1),[y*1.05,y*1.1],'-k', 'LineWidth',1);%left edge drop
plot([1;1]*xs(2),[y*1.05,y*1.1],'-k', 'LineWidth',1);%right edge drop
hold off
%--------------------------------------------------------------------------
% Helper function that Returns the largest single value of ydata in a given
% graphic handle. It returns the given value if it is not a graphics
% handle. 
function y=gety(h)
    %Returns the largest single value of ydata in a given graphic handle
    %h= figure,axes,line. Note that y=h if h is not a graphics
    if isgraphics(h) 
        switch(get(h,'type'))
            case {'line','hggroup','patch'},
                y=max(get(h,'ydata'));
                return;
            otherwise
                ys=[];
                hs=get(h,'children');
                for n=1:length(hs)
                    ys=[ys,gety(hs(n))];
                end
                y=max(ys(:));
        end
    else
        y=h;
    end
end
end