function download_toolboxes(toolbox_urls_file, download_dir)
% Downloads toolboxes from a list of URLs in a text file to a specified directory

% Throw an error if toolbox_urls_file is not defined
if ~exist('toolbox_urls_file', 'var') || isempty(toolbox_urls_file)
    toolbox_urls_file ='toolbox_urls.txt';
end

% Prompt the user to enter the download directory if it's not defined
if ~exist('download_dir', 'var') || isempty(download_dir)
    download_dir = input('Enter download directory: ', 's');
end

% Create the download directory if it doesn't already exist
if ~exist(download_dir, 'dir')
    mkdir(download_dir);
end

% Read the list of toolbox URLs from the text file
urls = importdata(toolbox_urls_file);

% Loop through each URL and download the corresponding toolbox
for i = 1:length(urls)
    url = urls{i};
    [~,filename,ext] = fileparts(url);
    download_path = fullfile(download_dir, filename);
    if strcmpi(ext, '.git')
        % Clone the Git repository if the directory doesn't exist
        if ~exist(download_path, 'dir')
            fprintf('Cloning %s to %s\n', url, download_path);
            command = sprintf('git clone %s %s', url, download_path);
            [status, cmdout] = system(command);
            if status ~= 0
                error('Error cloning %s: %s', url, cmdout);
            end
        else
            fprintf('%s already exists, skipping\n', download_path);
        end
    else
        % Download the file using websave
        fprintf('Downloading %s to %s\n', url, download_path);
        websave(download_path, url);
    end
end





end

