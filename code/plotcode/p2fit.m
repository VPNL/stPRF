function [yfit]= p2fit(x,y,xfit)
    [p,S] = polyfit(x,y,2); 
    yfit = polyval(p,xfit,S);

end