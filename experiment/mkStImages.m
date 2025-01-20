
function selected_images =  mkStImages(image_dir,nstim)

image_list = imageDatastore(image_dir,'IncludeSubfolders',true); % get db of images

fileNames = image_list.Files;
randomlySelectedFiles = {};  % Initialize an empty cell array to store selected files

for i = 1:length(nstim)
    pattern = sprintf('loc_%d_img', i);
    containsPattern = contains(fileNames, pattern);
    filteredFiles = fileNames(containsPattern);
    randomIndices = randperm(numel(filteredFiles), nstim(i));
    randomlySelectedFiles = [randomlySelectedFiles; filteredFiles(randomIndices)];
end
    selected_images = randomlySelectedFiles;
end

