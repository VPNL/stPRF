function setConstants(new_output_dir,new_simulator_dir,new_sessionDate)
% Read the contents of the original .m file
file_name = './template/c_template.m';  % Replace with the name of your original .m file
fid = fopen(file_name, 'r');
file_contents = fread(fid, '*char')';
fclose(fid);

output_dir = 'p1';
simulator_dir = 'p2';
sessionDate = 'p3';

% Modify the directory paths in the file contents
new_contents = strrep(file_contents, output_dir, new_output_dir);
new_contents = strrep(new_contents, simulator_dir, new_simulator_dir);
new_contents = strrep(new_contents, sessionDate, new_sessionDate);

% Write the modified file contents to a new .m file
new_file_name = 'code/Constants.m';  % Replace with the name of the new .m file
fid = fopen(new_file_name, 'w');
fwrite(fid, new_contents, 'char');
fclose(fid);

fprintf('global variable is set \n')
end