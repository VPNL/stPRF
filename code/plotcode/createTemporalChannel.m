function f = createTemporalChannel(params,modelNumber,windowSize)


if isempty(params)
    params.tau1 = 0.05; params.weight = 0; params.tau2 = 0.1;
    params.n = 2; params.sigma = 0.1; params.scale = 1;
    params.tau_s = 5;
end

params.fs = 100;
t    = 0:(1/params.fs):windowSize;

if modelNumber ==1
    f.temporal(:,1) = ones(1, length(t));
elseif modelNumber ==2
    onRange = 1:50;
    stim = zeros(1, length(t));
    stim(1, onRange) = 1;
    result = DNmodel(params, stim', 0);
    f.temporal(:,1) = result;
    
elseif modelNumber == 3
    temporal = zeros(length(t),3);
    
    n1 = 9; n2=10; kappa=1.33; fs =100; f.temporal=[];
    tau_s = params.tau_s;
    irfSustained = tch_irfs('S', tau_s, n1, n2, kappa, fs);
    irfTransient = tch_irfs('T', tau_s, n1, n2, kappa, fs);
    % f.temporal = zeros(max([length(irfSustained),length(irfTransient)]),3);
    % Normalize such that area under the curve sums to 1
    temporal(1:length(irfSustained),1) = nearZerotoZero(normSum(irfSustained));
    % For transient nrf, we do this separately for positive and negative parts:
    % First find indices
    pos_idx = irfTransient>=0;
    neg_idx = irfTransient<0;
    scf = 0.5;% 0.5; % sum of area under each pos/neg curve
    % Get positive part and normalize sum
    irfT_pos = irfTransient(pos_idx);
    irfT_pos = scf*normSum(irfT_pos);
    % Get negative part and normalize sum
    irfT_neg = abs(irfTransient(neg_idx));
    irfT_neg = -scf*normSum(irfT_neg);
    % Combine positive and negative parts
    nrfT2 = NaN(size(irfTransient));
    nrfT2(pos_idx) = irfT_pos;
    nrfT2(neg_idx) = irfT_neg;
    temporal(1:length(nrfT2),2) = nearZerotoZero(nrfT2);
    temporal(1:length(nrfT2),3) = nearZerotoZero(-nrfT2);
    f.temporal = temporal;
end