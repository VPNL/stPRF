
function  varexp = calculateR2(raw,prediction)

res = raw' - prediction';
res_var = sum(res .^ 2) ./ sum((raw' - mean(raw')) .^ 2);
varexp = 1 - res_var;

end
