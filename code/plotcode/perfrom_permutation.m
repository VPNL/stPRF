function p_value = perfrom_permutation(param_to_test,roiList)
%% permutation test
rng("default")

% param_to_test=window;
n_voxel = 100;
n_boot  = 100;
n_permute = 10000;
all_sample =  zeros(n_boot,length(roiList));
f_null = zeros(1,n_permute);
voxel_sample = zeros(n_voxel,length(roiList));
f_dist = zeros(1,n_permute);
for p = 1:n_permute
    for i = 1:n_boot        
        for r = 1:length(roiList)
            w = (cell2mat(param_to_test(:,r)));
            voxel_sample(:,r) = (datasample(w,n_voxel,'replace',true));
        end

        if p == 1 % get real data
            all_sample(i,:) = nanmedian((voxel_sample));
        else % get shuffled data
            all_sample(i,:) = nanmedian(shuffle_matrix(voxel_sample));
        end
    end

%%
stats_table1 = array2table(all_sample, 'VariableNames', roiList);
roi = table(roiList', 'VariableNames',{'roi'});
rm = fitrm(stats_table1, 'V1-IPS~1', 'WithinDesign',roi);
ranovaResults = ranova(rm);
f_dist(p) = ranovaResults.F(1);

% it does not really matter which stats you use... 
% you can also do lmm 
% rois = repmat(roiList,size(all_sample,1),1);
% T = table(all_sample(:), rois(:),'VariableNames',["value1","roi"]);
% lm = fitlm(T,'value1~roi');
% lm_anova = anova(lm);
% f_null(p) = lm_anova.F(1);

% check progress
% if mod(p,100) == 0
%     sprintf(".")
% elseif p == n_permute
%     sprintf("stats done\n")
% end
% p_value = sum(f_dist >= f_dist(1)) / size(f_null,2);

end

