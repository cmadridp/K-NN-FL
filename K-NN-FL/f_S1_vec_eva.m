function result = f_S1_vec_eva(x)
    % Determine the size of the input matrix x
    d = size(x, 2); % Number of columns in x, representing the dimensionality of the data
    n = size(x, 1); % Number of rows in x, representing the number of data points

    % Initialize a result vector of zeros with the same number of rows as x
    result = zeros(n, 1);

    % Create a row vector of ones with the same number of columns as x
    ones_1 = ones(1, d);  % Note: This vector is created but not used in the function

    % Loop over each row (data point) in x
    for i = 1:n
        % Apply the function f to each row of x and store the result
        result(i) = f(x(i, :));
    end
end
