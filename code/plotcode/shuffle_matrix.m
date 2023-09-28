function shuffled_matrix = shuffle_matrix(original_matrix)
    % Get the dimensions of the matrix
    [num_rows, num_columns] = size(original_matrix);

     shuffled_matrix = original_matrix(:,randperm(num_columns));
end
