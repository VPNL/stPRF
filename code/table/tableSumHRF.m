function DT= tableSumHRF(hrfTable)
DT=[];

rois = unique(hrfTable.name);
for s = unique(hrfTable.subj)'
    for r = 1:length(rois)
        params.subj=s;
        params.name=rois(r);
        newDT = tableFilter(hrfTable,params);
        
        tmpDT = [];
        tmpDT.subj = params.subj;
        tmpDT.name = params.name;
        tmpDT.coords = {newDT.coords};
        tmpDT.params = nanmean(newDT.params);
        tmpDT.hrfs =  nanmean(newDT.hrfs);
        tmpDT =  struct2table(tmpDT);
        
   
       DT = [DT ; tmpDT];
        
    end
end
end