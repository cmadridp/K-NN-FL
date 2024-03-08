function result = f(x)
    % Determine the number of columns in x
    d = size(x, 2);

    % Create a row vector of ones with the same number of columns as x
    ones_1 = ones(1, d);

    % Calculate the squared Euclidean norm (L2 norm) of (x - 1/3*ones_1)
    % This computes the sum of squares of the elements in x minus 1/3
    norm1 = sum((x - 1/3*ones_1).^2);

    % Calculate the squared Euclidean norm (L2 norm) of (2*x - 5/3*ones_1)
    % This computes the sum of squares of the elements in 2*x minus 5/3
    norm2 = sum((2*x - 5/3*ones_1).^2);

    % Compare the two norms and assign the result
    % If norm1 is smaller than norm2, result is set to 1, else result is -1
    if norm1 < norm2
        result = 1;
    else
        result = -1;
    end
end
