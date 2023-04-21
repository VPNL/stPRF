function [fn,saveName,tbl1] = st_plot_varexp_violin(DT, params,roiList)



% initial filter 
df = tableFilter(DT,params);
df = st_tablegiveROInumber(df, roiList);

[G,~] = findgroups(df.roiNumber,df.tmodel,df.subjID);
df.G = G;
eachGroup=[]; N = []; newT=[];
% edges = linspace(0,1,nbin);
for i = 1:length(unique(G))
    params.G = i;
    eachGroup = tableFilter(df,params);
    newT(i).roiNumber = double(eachGroup.roiNumber(1));
    newT(i).tmodel = eachGroup.modelNumber(1);
    
    subjID = str2num(regexp(eachGroup.subjID(1),'\d*','Match'));

    newT(i).subjID = subjID;
    newT(i).varExp = double( mean(eachGroup.cv_varexp));
    newT(i).all = double( (eachGroup.cv_varexp));

end
newT = struct2table(newT);
%%
[newT.G,~] = findgroups(newT.roiNumber,newT.tmodel);
mycolor = getmycolors(1);
g_color = repmat([mycolor; 0 0 0],length(unique(newT.G)),1);
for x = 1:length(roiList)
    idx(x) = (4*x)-1;
    newT.G(newT.G>idx(x)) = newT.G(newT.G>idx(x))+1;
end

%
% lme = fitlme(newT,'varExp ~ roiNumber + (1|subjID)');
%%
params = [];
st = []; g=[]; maxval=[];
for er = 1:length(roiList) 
    for t = 1:3
        params.roiNumber = er;
        params.tmodel = t;
        df = tableFilter(newT,params);
        g{er,t} = df.varExp;
%         g{er,t} = cell2mat(df.all);
    end
    
%     p(er) = signrank(g{er,2},g{er,3})
    p1 =fppermutationtestmp(g{er,1},g{er,2},1);
    p2 =fppermutationtestmp(g{er,1},g{er,3},1);
    p3 =fppermutationtestmp(g{er,2},g{er,3},1);
    
    st= [st; p1(3) p2(3) p3(3)];
    maxval(er) = max(max(cell2mat(g(er,:))));
end

[st]=fdr_bh(st);
sl = siglevel(st);
tbl1 = array2table([st sl],"VariableNames", ...
    ["c1","c2","c3","sig1","sig2","sig3"]);



%%
fn = tch_fig('varexp',[12.5413    6.0325   22.5425    8.2021]);
saveName = 'violin-var_exp';

vs = violinplot([newT.varExp; ones(size(idx))'*-.01], ...
    [newT.G; idx'+1], ...
    'ViolinColor',g_color, ...
    'ViolinAlpha', 1, ...
    'MedianColor',[0 0 0], ...
    'ShowData',false, ...
    'ShowBox', false, ...
    'ShowWhiskers',false);

ylabel('variance explained (R^2)')
xticks(idx-1); yticks(0:0.2:1); ylim([0 0.8]);
xticklabels(roiList);
tch_set_axes;

%%
for esl = 1 :size(sl,1)
    pos = 1+(esl-1)*4;

    if sl(esl,1) > 0                   
        sigline([pos,pos+1],repmat('*',1,sl(esl,1)),maxval(esl)); hold on
    end

    if sl(esl,2) > 0
        sigline([pos,pos+2],repmat('*',1,sl(esl,1)),maxval(esl)+maxval(esl)*0.15)
    end
    
    if sl(esl,3) > 0
        sigline([pos+1,pos+2],repmat('*',1,sl(esl,1)),maxval(esl))
    end
    
end



%%


    function  [pValues] = fppermutationtestmp(G1Data,G2Data,simulate,permutations)
        % Produces p-values from a Fisher-Pitman permutation test of matched-pairs.
        % Input arguments
        %
        % G1Data - A column vector of the first group of data
        % G2Data - A column vector of the second group of data with paired observations
        % to group 1.
        % simulate - set to 1 to simulate the full permutation test
        % permutations - the number of permutations to run if simulation is chosen
        % (default 50000)
        % Check that the Groups are the same size
        if any(size(G1Data) ~= size(G2Data)) || size(G1Data,2)~=1;
            error('Data vectors must be column vectors of the same size')
        end
        %Revert to default number of permutations if unspecified
        if nargin == 3
            permutations = 50000;
        end
        %Initialise variables
        n = size(G1Data,1); %number of pairs
        DifferenceVector = G1Data-G2Data;
        if simulate == 1
            SignMatrix = 1 - 2*round(rand(permutations,n));
        else
            permutations = 2^n;
            
            if permutations > 1000000
                check = input(['Running with ',num2str(permutations) ,' permutations. Continue? [Y/N]: '],'s');
                if check =='N'
                    pValues = [];
                    return
                end
            end
            
            D = (0:permutations-1)';
            SignMatrix = 1 - 2*rem(floor(D*pow2(-(n-1):0)),2);
        end
        % Compute the statistic for all the permutations:
        StatDist = SignMatrix*DifferenceVector;
        % Critical value for given allocation
        critvalue = sum(DifferenceVector);
        pValues(1) = sum(StatDist>=critvalue,1);
        pValues(2) = sum(StatDist<=critvalue,1);
        pValues(3) = 2*min(pValues);
        pValues = pValues/permutations;
    end

    function [adj_p]=fdr_bh(pvals)
        
        q=.05;
        method='pdep';
        
        s=size(pvals);
        if (length(s)>2) || s(1)>1,
            [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
        else
            %p-values are already a row vector
            [p_sorted, sort_ids]=sort(pvals);
        end
        [~, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
        m=length(p_sorted); %number of tests
        if strcmpi(method,'pdep'),
            %BH procedure for independence or positive dependence
            threshvals=(1:m)*q/m;
            wtd_p=m*p_sorted./(1:m);
            
        elseif strcmpi(method,'dep')
            %BH procedure for any dependency structure
            denom=m*sum(1./(1:m));
            threshvals=(1:m)*q/denom;
            wtd_p=denom*p_sorted./[1:m];
            %Note, it can produce adjusted p-values greater than 1!
            %compute adjusted p-values
        else
            error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
        end
        %compute adjusted p-values; This can be a bit computationally intensive
        adj_p=zeros(1,m)*NaN;
        [wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
        nextfill = 1;
        for k = 1 : m
            if wtd_p_sindex(k)>=nextfill
                adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
                nextfill = wtd_p_sindex(k)+1;
                if nextfill>m
                    break;
                end;
            end;
        end;
        adj_p=reshape(adj_p(unsort_ids),s);
        
    end

    function  assigned_values = siglevel(A)
        tv = [0.05, 0.01, 0.001];
        assigned_values = zeros(size(A));
        for ts = 1:length(tv)
            k  = A < tv(ts);
            assigned_values = assigned_values + k;
        end
    end
end