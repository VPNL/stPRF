function model = updateModel(df,model,fold,params)

% fields = fieldnames(model{1});
% tm =  model{1};
% for ef = 1:length(fields)
%     fields{ef}
% end

tmp = df(fold);

model{1}.x0 = tmp.x0;
model{1}.y0 = tmp.y0;
model{1}.sigma = tmp.sigma;
model{1}.exponent = tmp.exponent;
model{1}.beta = tmp.beta;
model{1}.varexpfitprf = tmp.varexp;


zArray = zeros(size(tmp.x0));
    switch lower(params.analysis.temporalModel)
        case {'glm','1ch-glm'}
        case {'dn','1ch-dcts'}                       
            model{1}.tau1       = zArray; 
            model{1}.weight     = zArray;
            model{1}.tau2       = zArray;
            model{1}.nn         = zArray;
            model{1}.delay      = zArray;
        case {'3ch','3ch-stln'}
            model{1}.tau_s      = zArray;
            model{1}.tau_t      = zArray;

    end

% 
% model{1}.varexp = zeros(size(model{1}.rawrss));
% model{1}.varexp(find(model{1}.rawrss)) = tmp.varexp;

 
end