function ari = adjusted_rand_index(comm1, comm2)
    % Create contingency table
    cont = accumarray([comm1(:), comm2(:)], 1);
    
    % Sum over rows and columns
    sum_rows = sum(cont, 2);
    sum_cols = sum(cont, 1);
    
    % Calculate combinations (n choose 2)
    n = sum(cont(:));
    total_pairs = n * (n - 1) / 2;
    
    % Compute a, b, c, d
    a = sum(cont(:) .* (cont(:) - 1)) / 2;
    sum_row_pairs = sum(sum_rows .* (sum_rows - 1)) / 2;
    sum_col_pairs = sum(sum_cols .* (sum_cols - 1)) / 2;
    
    expected_a = sum_row_pairs * sum_col_pairs / total_pairs;
    max_a = (sum_row_pairs + sum_col_pairs) / 2;
    
    % Adjusted Rand Index
    ari = (a - expected_a) / (max_a - expected_a);
end