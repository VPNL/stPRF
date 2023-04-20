function randjson(jsonDir,templateJson,nVoxels,expNames,stimFiles,StimSize,temporalSampleRate,RF,tParam1,tParam2,noiselevel)
% create json files with spatiotmpoeral parameters being specified

if length(expNames) ~= length(stimFiles)
    error("expName and stimFile length must match")
end


    J = loadjson(templateJson);

for jj = 1:length(expNames)
    
    J = loadjson(templateJson);

    savename{1} =  convertStringsToChars(expNames(jj));

%     [~,jsonname]=fileparts(templateJson{jj});
%     J = loadjson(templateJson{jj});
    
    if iscell(J)
        J=J{:};
    end
    
    jsonData = J;
    
    % check repeats
    jsonData.repeats = 1;
    
    
    
    % check format
    vs = {'Temporal','HRF','RF','Noise','Stimulus','Temporal'};
    for ev = 1:length(vs)
        if ~iscell(jsonData.(vs{ev}))
            
            % Get the size of the struct array
            [~,n] = size(jsonData.(vs{ev}));

            % Create a cell array to hold the extracted structs
            tempCell = cell(1, n);

            % Loop through the struct array and extract each struct into a 1x1 cell array
            for i = 1:n
                tempCell{i} = jsonData.(vs{ev})(i);
            end

            % Assign the cell array to jsonData.Temporal
            jsonData.(vs{ev}) = tempCell;
        end
    end

    
    % set noise level and save name
    if strcmp(noiselevel,'noise1') ||  strcmp(noiselevel,'noise2') ||  strcmp(noiselevel,'low')
        jsonData.Noise{1}.voxel = "low"; % low noise added
        jsonData.Noise{1}.seed = "random";
        savename{3} = sprintf('%s.json',noiselevel);
    elseif strcmp(noiselevel,'noise0') ||  strcmp(noiselevel,'none')
        jsonData.Noise{1}.voxel = "low"; % Note that this is a dummy. actually no-noise is added
        jsonData.Noise{1}.seed = "none";
        savename{3} = sprintf('%s.json',noiselevel);        
    elseif  isnumeric(noiselevel)
        jsonData.Noise{1}.voxel = "SNR"; % Note that this is a dummy. Noise level will be matched to the value
        jsonData.Noise{1}.seed = "random";
        savename{3} = sprintf('SNR%d.json',noiselevel);
    end
    
    % stimulus:
    jsonData.Stimulus{1}.fieldofviewHorz = StimSize;
    jsonData.Stimulus{1}.fieldofviewVert = StimSize;
    jsonData.Stimulus{1}.expName = convertStringsToChars(expNames(jj));
    jsonData.Stimulus{1}.myload  = convertStringsToChars(stimFiles(jj));
    if ~isfile(jsonData.Stimulus{1}.myload)
        jsonData.Stimulus{1}.myload
        error('stimulus file does not exist')
    end
    %% assign params
    for vv = 1:nVoxels
        jsonData.RF{1}.Centerx0    =  RF(vv,1);
        jsonData.RF{1}.Centery0    =  RF(vv,2);
        jsonData.RF{1}.sigmaMajor  =  RF(vv,3);
        jsonData.RF{1}.sigmaMinor  =  RF(vv,3);
        
        for tt = 1:length(jsonData.Temporal)
            jsonData.Temporal{tt}.fs = temporalSampleRate;
            switch jsonData.Temporal{tt}.temporalModel
                case {'3ch-stLN', 'CST'}
                    jsonData.Temporal{tt}.tParams.values = tParam1(vv,:);
                case {'1ch-dcts', 'DN-ST'}
                    jsonData.Temporal{tt}.tParams.values = tParam2(vv,:);
            end
        end

        % update Voxelnumber as savename
        savename{2} = ['voxels' num2str(sprintf('%03d', vv))];
        savethisname = strjoin(savename,'_');
        
        % clearn up and change
        jsonString = jsonencode(jsonData);
        jsonString = strrep(jsonString, ',', sprintf(',\n'));
        jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
        jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));
        jsonString = jsonString;
        fid = fopen(fullfile(jsonDir,savethisname), 'w');
        fwrite(fid, jsonString,'char');fclose(fid);
        %
        
    end
    
end
end