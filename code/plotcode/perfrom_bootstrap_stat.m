function ranovaResults = perfrom_bootstrap_stat(param_to_test,roiList)
%% boostrap stat test
rng("default")

% param_to_test=window;
n_voxel = 1000;
n_boot  = 100;
all_sample =  zeros(n_boot,length(roiList));
voxel_sample = zeros(n_voxel,length(roiList));
for i = 1:n_boot
    for r = 1:length(roiList)
         w = (cell2mat(param_to_test(:,r)));
        voxel_sample(:,r) = (datasample(w, n_voxel,'replace',true));
    end
    all_sample(i,:) = nanmedian((voxel_sample));
end



%%
stats_table1 = array2table(all_sample, 'VariableNames', roiList);
roi = table(roiList', 'VariableNames',{'roi'});
rm = fitrm(stats_table1, 'V1-IPS~1', 'WithinDesign',roi);
ranovaResults = ranova(rm);
% f_dist(p) = ranovaResults.F(1);

% it does not really matter which stats you use... 
% you can also do lmm 
% rois = repmat(roiList,size(all_sample,1),1);
% T = table(all_sample(:), rois(:),'VariableNames',["value1","roi"]);
% lm = fitlm(T,'value1~roi');
% lm_anova = anova(lm);

end



