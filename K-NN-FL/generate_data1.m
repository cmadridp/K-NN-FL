function [X_s, y_i_js, X_s_training, y_i_js_training, Weights_training, X_s_test, y_i_js_test, Weights_test, Weights] = generate_data1(n, d, mult)
    % Initialize multiplier-based parameters
    m1 = 10 * (mult + 2);
    m2 = 10 * (mult);
    m3 = 10 * (mult + 3);
    m4 = 10 * (mult + 1);

    % Define the n_ts
    ms = NaN(n, 1);  % Initialize an array to store 'ms' values
    % Assign different 'mi' values to quarters of the 'ms' array
    for i = 1:(n / 4 + 1)
        ms(i) = m1;
    end
    for i = (n / 4 + 1):(n / 2 + 1)
        ms(i) = m2;
    end
    for i = (n / 2 + 1):(3 * n / 4 + 1)
        ms(i) = m3;
    end
    for i = (3 * n / 4 + 1):n
        ms(i) = m4;
    end
    nm = sum(ms);  % Total count across all 'ms' values

    % Creating the X_s matrix
    X_s = NaN(nm, d);  % Initialize matrix for feature vectors
    % Generate random feature vectors
    for i = 1:nm
        for j = 1:d
            X_s(i, j) = rand(1);  % Generate a uniform random number between 0 and 1
        end
    end

    % Creating the delta_i / b_t_is matrix
    b_i_ts = NaN(25, nm);  % Initialize matrix for b_i_ts values
    % Generate random b_i_ts values
    for i = 1:nm
        for t = 1:25
            b_i_ts(t, i) = normrnd(0, 1);  % Generate a normally distributed random number with mean 0 and standard deviation 1
        end
    end

    % Creating the epsilon_s vector
    epsilon_s = NaN(nm, 1);  % Initialize vector for epsilon_s values
    % Generate random epsilon_s values
    for i = 1:nm
        epsilon_s(i, 1) = randn(1, 1);  % Generate a normally distributed random number with mean 0 and standard deviation 1
    end

    % Creating the y_i_js and Weights vectors
    y_i_js = NaN(nm, 1);  % Initialize vector for y_i_js values
    Weights = NaN(nm, 1);  % Initialize vector for Weights
    aux = 0;  % Initialize auxiliary variable
    % Generate y_i_js and Weights values
    for i = 1:n
        aux = ms(i) + aux;
        for j = 1:ms(i)
            Weights(aux - ms(i) + j) = 1 / (n * ms(i));
            y_i_js(aux - ms(i) + j) = f(X_s(aux - ms(i) + j, :)) + xi_i(X_s(aux - ms(i) + j, :), i, d, b_i_ts) + epsilon_s(aux - ms(i) + j);
        end
    end

    % Create training and test indices
    train_indices = datasample(1:length(y_i_js), round(length(y_i_js) * 0.75), 'Replace', false);

    % Subset the data into training and test sets
    X_s_training = X_s(train_indices, :);
    y_i_js_training = y_i_js(train_indices);
    Weights_training = Weights(train_indices);
    X_s_test = X_s(setdiff(1:length(y_i_js), train_indices), :);
    y_i_js_test = y_i_js(setdiff(1:length(y_i_js), train_indices));
    Weights_test = Weights(setdiff(1:length(y_i_js), train_indices));
end
