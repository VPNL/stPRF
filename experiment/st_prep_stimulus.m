function st_prep_stimulus()
%% experimental setting
% stimulus radius = 9
% trial Dur = 5, blank Dur =10
% Number of Direction 4, number of steps 9
% Stimulus [4 * 9 * 5] = 180 sec
% Blank (10) + Stimulus (4 * 9 * 5 = 180) + 2 Blanks (10*2)+ Blank (10)
% Total duration = 220 sec => 3.6 min
% load predefined shuffled sequence for each run

% save run Random sequence & parfiles

%% SET PATHS
cd  '/Users/insubkim/Documents/MATLAB/exp/stExp'
addpath(genpath(    '/Users/insubkim/Documents/MATLAB/Psychtoolbox-3'))
addpath(genpath(    '/Users/insubkim/Documents/MATLAB/exp/stExp'))
% addpath(genpath('/Users/insubkim/Desktop/scanning_CNI/toonExp'))
%% 
%% 
% addpath(genpath('/Users/insubkim/Documents/vistasoft'))

%% constants

%subject
[a b]=GetKeyboardIndices;

deviceNumber= max(a);
subject.name = 'chop';

subject.scanner = 0;
% stimulus
params.radius = 12;
params.doShuffle = 0;

% temporal
params.trial_dur = 5;
params.blank_dur = 10;
params.tr = 1;

% scanner
triggerKey = 't';
comTrigScan = true;
% comTrigScan = false;

%save
saveStim = 0;

%eyelink
params.eyelink = 0;


%% COLLECT SESSION INFORMATION
% initialize subject data structure
subject.name = deblank(subject.name);
subject.date = date;


while ~ismember(params.eyelink,[0 1])
    params.eyelink = input('Would you like to eyetrack? (Y=1,N=0): ');
end

if params.eyelink
    fprintf(['\nwfret paused to make sure that you have turned off the' ...
             ' wifi, press any key to continue\n\n'])
    pause;
end



%%
path.baseDir = pwd; addpath(genpath(path.baseDir));
path.dataDir = fullfile(path.baseDir,'data',subject.name);

if ~exist(path.dataDir) &&  ~strcmp(subject.name,'tes'), mkdir(path.dataDir); end
runN = getAllFiles(path.dataDir,'*.txt',1);
runN=length(runN)+1;



if runN ==1 && exist([path.dataDir '/ExpDesign.mat'], 'file')  ~= 2
    ExpDesign =  mkExpDesign(36,9);
    save([path.dataDir '/ExpDesign.mat'], 'ExpDesign');
    disp('######## Created and Saved ExpDesign');
else
   load([path.dataDir '/ExpDesign.mat']);
   disp('######## loading ExpDesign');
end


% subject.temptype = ExpDesign{runN,1};
% if strcmp(ExpDesign{runN,2},'seq')
%     subject.doshuffle = 0;
% else
%     subject.doshuffle = 1;
% end

subject.nRun =runN;
subject.saveName = fullfile(path.dataDir,sprintf('stPRF_Run%d.txt',runN));
subject.ExpDesign =ExpDesign;

 
ExpDesign = ExpDesign(runN,:);

% subject.temptype = 'c';
% subject.doshuffle = 0;
% 

%% Experiment params

% load display
params.display =loadDisplayParams;


%% override - HARVARD 
params.display.backColorIndex = 176; %(0.69)
ppd = 34.4284;




temporal.trial_dur = params.trial_dur;
temporal.fs = 1/params.display.frameRate;
viewTime = temporal.fs;

params.temporal.motionSteps = 30;

params.stimSize = params.radius;
params.temporal.frequency =1;
params.ncycles = 1;
params.period = 36;
params.ringDeg          = params.radius./4;			% Ring radius/width (deg)
params.prescanDuration = 0;
params.framePeriod  = params.tr;
params.numImages    = round(params.period/params.framePeriod)*2;
params.type = 'bar';

% various time measurements:
duration.stimframe          = 1./params.temporal.frequency./(params.temporal.motionSteps);
duration.cycle.seconds      = params.period;


params.ringWidth        = params.ringDeg;
ringWidth       = params.ringWidth;
halfNumImages   = params.numImages./2;
numMotSteps     = 63 ;%params.temporal.motionSteps;
outerRad=params.radius;

numSubRings     = 1;
bk = params.display.backColorIndex;

% made outerRad to be 1/3 because matrix becomes to big
% I will re-scale the image it later on
m = round(2 * ppd* outerRad);
n = round(2 * ppd* outerRad);
[x,y]=meshgrid(linspace(-outerRad,outerRad,n),linspace(outerRad,-outerRad,m));

% r = eccentricity;
r = sqrt (x.^2  + y.^2);

% loop over different orientations and make checkerboard
% first define which orientations
orientations = (0:45:360)./360*(2*pi); % degrees -> rad
orientations = orientations([1 6 3 8 5 2 7 4]);
remake_xy    = zeros(1,params.numImages)-1;
remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;
original_x   = x;
original_y   = y;
% step size of the bar
step_nx      = duration.cycle.seconds./params.tr/4;
step_x       = (2*outerRad) ./ step_nx;
step_startx  = (step_nx-1)./2.*-step_x - (ringWidth./2);
%[0:step_nx-1].*step_x+step_startx+ringWidth./2
fprintf('[%s]:stepsize: %f degrees.\n',mfilename,step_x);

% if we create colored bars we want to make the edges soft.
softmask = ones(m);

% Loop that creates the final images
fprintf('[%s]:Creating %d images:',mfilename,halfNumImages);
images=zeros(m,n,3,halfNumImages*numMotSteps,'uint8');
%%
picked = randperm(6,1);
for picked = 1
cacheImg =[path.baseDir '/img_' num2str(picked) '.mat'];
if exist(cacheImg, 'file')  ~= 2
    for imgNum=1:halfNumImages
        
        if remake_xy(imgNum) >=0,
            x = original_x .* cos(remake_xy(imgNum)) - original_y .* sin(remake_xy(imgNum));
            y = original_x .* sin(remake_xy(imgNum)) + original_y .* cos(remake_xy(imgNum));
            % Calculate checkerboard.
            % Wedges alternating between -1 and 1 within stimulus window.
            % The computational contortions are to avoid sign=0 for sin zero-crossings
            wedges    = sign(round((cos((x+step_startx)*numSubRings*(2*pi/ringWidth)))./2+.5).*2-1);
            posWedges = find(wedges== 1);
            negWedges = find(wedges==-1);
            rings     = zeros(size(wedges));
            
            checks    = zeros(size(rings,1),size(rings,2),params.temporal.motionSteps);
            
            % reset starting point
            loX = step_startx - step_x;
        end;
        
        avail_imgs = 63;
        potential_imgs = randperm(avail_imgs,numMotSteps);
        for ii=1:numMotSteps,
            
            % load image
            % update the image to be the same size
            filename = ['pic' mat2str(potential_imgs(ii)) '.jpg'];
            newfile = imread(filename);
            cropped_img = imresize(newfile,[m n]);

%             newfile = imresize(newfile,[m/2 n/2]);
%             cropped_img = imresize(newfile,0.2);
%             cropped_img = imresize(cropped_img,5);

%             cropped_img = newfile(1:m, 1:n,:);

            toons(:,:,:,ii) = cropped_img;
            
        end;
        
        loX   = loX + step_x;
        hiX   = loX + ringWidth;
        
        window = ( (x>=loX & x<=hiX) & r<outerRad);
        
        
        tmpvar = zeros(m,n);
        tmpvar(window) = 1;
        tmpvar = repmat(tmpvar, [1 1 3]);
        tmpvar = repmat(tmpvar,[1 1 1 numMotSteps]);
        window = tmpvar == 1;
        img         = bk*ones(size(toons));
        img(window) = toons(window);
        images(:,:,:,(imgNum-1).*numMotSteps+1:imgNum.*numMotSteps) = uint8(img); %#ok<*BDSCA>
        
        
        
        fprintf('.');drawnow;
        
    end
    save(cacheImg,...
        'images', '-v7.3');

else
    load(cacheImg)
end

fprintf('Done.\n');
end

% %% cut it and save it
% for picked = 1
%     cacheImg =['./img_' num2str(picked) '.mat'];
%     load(cacheImg)
%     mkdir(sprintf("./stPRF_img/imgset%d/",picked))
%     for di = 1:36
%         for  i = 1:numMotSteps
%             imwrite(images(:,:,:,i),sprintf("./stPRF_img/imgset%d/loc_%d_img_%d.png",picked,di,i));
%         end
%         images(:,:,:,1:numMotSteps) = [];
%     end
% end
% %% get back image size
% m = round(2 * ppd* outerRad);
% n = round(2 * ppd* outerRad);
% 
% images = imresize(images,[m n]);


%% actual sequence

[img_sequence, nstim]= mkTemporal(ExpDesign,temporal);

temp = 1:size(images,4);
samp = reshape(temp,numMotSteps,[]); 

imageNumber=[];
for q = 1:params.period % 36
    imageNumber{q} = samp(1:nstim{q},q);
end
imageNumber = cellfun(@transpose, imageNumber,'UniformOutput',false);


sequence = [];

for q= 1:length(imageNumber)
    for i = 1:length(imageNumber{q})
        curpuls = imageNumber{q};
        sequence = [sequence curpuls(i)*img_sequence{q}];
    end

end

% sequence [frames X spatial location] 
% right now [300 frames = 5 sec, trial dur X 36 spatial location (9 step X 4 dir)] 
nF = params.trial_dur*params.display.frameRate;
nbF =  params.blank_dur*params.display.frameRate;
sequence = reshape(sequence,nF,[]);


% add Blanks
insert = @(a, x, n) cat(2,  x(1:n), a, x(n+1:end));
insert2 = @(a, x, n) cat(2,  x(:,1:n), a, x(:,n+1:end));
nImages = cell2mat(nstim);

% add blank inbetween
blankInd = [9 20];
for i = 1:length(blankInd)
    sequence  = insert2(zeros(nF,nbF/nF), sequence, blankInd(i));
    nImages   = insert(zeros(1,2), nImages, blankInd(i));
end

%add blank in the start & end
nImages  = [0 0 nImages 0 0];
sequence = [zeros(nF,nbF/nF) sequence zeros(nF,nbF/nF)];
sequence = sequence(:);
rs = reshape(sequence,nF,[]);


%% save stim
if saveStim == 1
    blankImage = uint8(ones(101,101,3).*bk);
    rs = reshape(sequence,nF,[]);
    stim=[];
    for a = 1:size(rs,2)
%         stim_trial = uint8(zeros(size(images,1),size(images,2),3,nF));
        stim_trial = uint8(zeros(101,101,3,nF));

        trial = rs(:,a);
        for i = 1:length(trial)
            if trial(i) ==0
                stim_trial(:,:,:,i) = blankImage;
            else
                curImg = images(:,:,:,trial(i));
                curImg = imresize(curImg,[101 101]);
                stim_trial(:,:,:,i) = curImg;
            end
        end
%         stim_trial = imresize(stim_trial,[101 101]);
        stim{a} = stim_trial;
    end
    save([path.dataDir '/stim_' num2str(runN) '.mat'],'stim','sequence','-v7.3');
end
%% Save lastbit
save([path.dataDir '/info_' num2str(runN) '.mat'],'expInfo','subject','data','params','-v7.3');

end
    