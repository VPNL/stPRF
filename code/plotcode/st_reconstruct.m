function df = st_reconstruct(DT,center,computeST)

if isempty(DT)
   return; 
end

if notDefined('center') 
    center = 0;
end

if notDefined('computeST') 
    computeST = 0;
end

windowSize = 5;
df = tableArrayFormat(DT);

%%
df.sRF = cell(height(df),1);
df.tRF = cell(height(df),1);
sRF =cell(height(df),1);
tRF =cell(height(df),1);

% spatial RF
tic
for erf = 1: height(df)
    if erf ==1
        fprintf('creating sRF: \n')
    end
    if rem(erf,1000) == 0
        fprintf('%d / %d \n',erf, height(df))
    end
    
    rf.sigmaMajor= df(erf,:).sigma;
    rf.sigmaMinor = rf.sigmaMajor;
    rf.Centerx0= df(erf,:).x0;
    rf.Centery0= df(erf,:).y0;

    % centered
    if center == 1
        rf.Centerx0= 0;
        rf.Centery0= 0;
    end
    
    rf.Theta = 0;
    sRF(erf) = {createRF(rf)};
end

toc
%%
tic
for erf = 1:height(df)
    if erf ==1
        fprintf('creating tRF: \n')
    end
    if rem(erf,1000) == 0
        fprintf('%d / %d \n',erf, height(df))
    end
    tf = df(erf,:).temporal;
    
    tparams=[]; targetModel = 1;
    if df(erf,:).tmodel == "1ch_dcts"
        targetModel = 2;
        tparams.tau1   = tf(1);
        tparams.weight = tf(2);
        tparams.tau2   = tf(3);
        tparams.n      = tf(4);
        tparams.sigma  = tf(5);
        tparams.scale  = 1;
    elseif df(erf,:).tmodel == "3ch_stLN"
        targetModel = 3;
        tparams.tau_s   = tf(1);
    end
    
    f = createTemporalChannel(tparams,targetModel,windowSize);
    if df(erf,:).tmodel == "1ch_dcts"
        f.temporal = nearZerotoZero((normSum(f.temporal)));
    end
    tRF(erf)= {f.temporal};        
end
toc

%%
% st RF
if computeST == 1
    for erf = 1:height(df)
        each_sRF = sRF{erf};
        [ ~, ix ] = max(each_sRF(:));
        [ i1,i2] = ind2sub( size(each_sRF), ix );
        px(:,erf) = each_sRF(i1,:);
        py(:,erf) = each_sRF(:,i2);

        stRF_sus{erf} =  py(:,erf)*tRF{erf}(:,1)';

        if size(tRF{erf},2) > 1
            stRF_tran{erf} = py(:,erf)*tRF{erf}(:,2)';
        end
    end
end


%% 
df.sRF = sRF;
df.tRF = tRF;

if computeST == 1
    df.px = num2cell(px',2);
    df.py = num2cell(py',2);
    df.stRF_sus = stRF_sus';
    if size(tRF{erf},2) > 1
        df.stRF_tran = stRF_tran';
    end
end


end




