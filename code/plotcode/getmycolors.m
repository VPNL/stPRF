function mycolor = getmycolors(colorSet)
% myCOLORsHEX = {'#009ADE','#C77CFF','#FF1F5B','#FFC61E','#F28522','#A0B1BA'};
if nargin < 1 || isempty(colorSet); colorSet = 1; end


if colorSet == 1
    %     myCOLORsHEX = {'#009ADE','#C77CFF','#FF1F5B'};
    myCOLORsHEX = {'#0E4BEF','#5a005a','#FF0090'};
    myCOLORsHEX = {'#0043ce','#6929c4','#9f1853'};
    myCOLORsHEX = {'#0072c3','#8a3ffc','#d02670'};
    myCOLORsHEX = {'#1e88e5','#ffc107','#d81b60'};
    myCOLORsHEX = {'#1e88e5','#f1c21b','#d81b60'};
    %     myCOLORsHEX = {'#514AC3','#f1c21b','#d81b60'};
    
elseif colorSet == 2
    myCOLORsHEX = {'#7570b3','#1b9e77' '#d95f02'  '#e7298a'};
elseif colorSet == 3
%     myCOLORsHEX = {'#7570b3','#584ed5', '#4335ee', ...
%         '#D81B60','#1b9e77', ...
%         '#d95f02','#9f4501',...
%         '#f393c4','#e7298a'};
%     
    myCOLORsHEX = {'#1ee5df','#1e88e5', '#1e25e5', ...
        '#00e600','#105b45  ', ...
        '#D81B60','#971343',...
        '#FFC107','#997404'};

    %     
    myCOLORsHEX = {'#1ee5df','#1e88e5', '#000066', ...
        '#00e600','#105b45  ', ...
        '#d81bbf','#C41E3A	',...
        '#FFC107','#997404'};
    
elseif colorSet == 4 % stream color
    %   mycolor = [0 0 1 ; 0.1 0.5 0.3;1 0.5 0;1 0 1];
    mycolor = [0 0 1; ...
        0.1 0.7 0.1 ; ...
        0.9 0.9 0 ;  ...
        0.9 0 0];
    
elseif colorSet == 5 % stream color
    mycolor = [0 0 0.3 ; 0 0 .6 ;0 0 1; ...
        0 0.5 0; 0.1 0.7 0.1 ; ...
        0.6 0 0 ; 0.9 0 0; ...
        0.9 0.7 0; 0.9 0.9 0 ;];
    
    
end

if colorSet < 4
    for cc = 1:length(myCOLORsHEX)
        mycolor(cc,:) = hex2rgb(myCOLORsHEX{cc});
    end
end
end