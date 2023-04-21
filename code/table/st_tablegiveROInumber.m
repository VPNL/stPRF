function df = st_tablegiveROInumber(df, roiList)

% give number to ROIs
pat = ("rh"|"lh");
for r = 1:length(roiList)
    if contains(roiList{r},pat) == 0
        df.roiNumber(find(df.roiName2==roiList{r})) = r;
    elseif  contains(roiList{r},pat) > 0
        df.roiNumber(find(df.roiName==roiList{r})) = r;
    end
end
% remove ROIs that are not used
df(find(df.roiNumber ==0),:) = [];

end