function stimulus = createStimulusStruct(img,cmap,seq,destRect,seqtiming,fixseq)%stimulus = createStimulusStruct(img,cmap,seq,destRect)%%Creates an stimulus structure with the fields:%%stimulus.image      array of images%stimulus.cmap       array of cmaps%stimulus.seq        array of indices into images and cmaps%stimulus.seqtiming  array of times associated with stimulus.seq%stimulus.fixSeq     binary array associated with fixation of seq%%stimulus.imagePtr		empty field used by Screen.%stimulus.srcRect		empty field used by Screen.%stimulus.destRect		where you want the stim to go%if nargin < 6,    fixseq = [];end;if nargin < 5,    seqtiming = [];end;if(~iscell(img)), img = {img}; endstimulus.images = img;stimulus.cmap = cmap;stimulus.seq = seq;stimulus.seqtiming = seqtiming;stimulus.fixSeq = fixseq;stimulus.textures = [];stimulus.srcRect = [];stimulus.destRect = [];return;