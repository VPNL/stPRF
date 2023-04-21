function DT= tableRemoveHRFDuplicateCoords(hrfTable)
DT=[];

for s = unique(hrfTable.subj)'
    for r =  unique(hrfTable.name)'
        params.subj=s;
        params.name=r;
        newDT = tableFilter(hrfTable,params);
        
        tmpDT = [];
        tmpDT.subj = repmat(params.subj,length(cell2mat(newDT.coords)),1);
        tmpDT.name = repmat(params.name,length(cell2mat(newDT.coords)),1);
        tmpDT.coords = cell2mat(newDT.coords);
        tmpDT.params = cell2mat(newDT.params);
        tmpDT.hrfs = cell2mat(newDT.hrfs);
        tmpDT.hrfs = cell2mat(newDT.hrfs);

        tmpDT =  struct2table(tmpDT);
        
        [~,idx]=unique(tmpDT.coords,'rows','first'); 
        tmpDT =  tmpDT(idx,:);
   
       DT = [DT ; tmpDT];
        
    end
end
end