function df= defaultThresh(DT,doBeta)

thresh=[];

if notDefined('doBeta')
    doBeta = 1;
end

thresh.varexp = 0.1;
df = tableThreshold(DT,thresh,'lower'); 

thresh = [];
thresh.ecc = 12;  % max ecc
df = tableThreshold(df,thresh,'upper');

if doBeta ==  1
    df = removeBetaOutlier(df);
end

end