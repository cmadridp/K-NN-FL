%% generate data
tic  % Start timer to measure execution time

result_rep=zeros(1,1);  % Initialize a matrix to store results of each repetition
for rep=1:1  % Loop for repetitions, currently set to run only once

    % Initialize variables for data generation
    n = 1000;  % Number of time points
    d = 2;     % Dimension of covariates
    mult = 0.5;  % Generator of m_i's (number of measurements per time)

    % Generate data using a custom function 'generate_data1'
    [X_s,y_i_js,X_s_training, y_i_js_training, Weights_training, X_s_test, y_i_js_test, Weights_test,Weights]=generate_data1(n,d,mult);

    %% Model Configuration and Cross-Validation Setup
    % Define a set of lambda values for regularization
    lambda=[0.00001,0.0001,0.01,0.1,1,10,100].';
    % Sort lambda values in ascending order
    new_lambda = sort(lambda);

    % Parameters for Cross-Validation
    folds_number = 5;  % Number of folds in cross-validation
    number_of_neighb = 5;  % Number of neighbors for kNN

    %% Cross-Validation to Find Optimal Lambda
    % Find the optimal lambda value using cross-validation
    lambda_opt=cv(X_s_training,ms_training,Weights_training,y_i_js_training,folds_number,number_of_neighb ,lambda,0.00001,5,10,2,2);
    %% Model Training with Optimal Lambda
    % Train the model using ADMM algorithm with k-NN-FL
    theta_s=admm_knnfl_varying_rho(X_s_training,y_i_js_training,Weights_training,lambda_opt, 5,0.00001, 5,10,2,2);

    %% Prediction on Test Data
    y_predicted = NaN(size(y_i_js_test));  % Preallocate matrix for predictions
    % Loop over each test sample
    for i = 1:size(y_predicted, 1)
        % Calculate squared Euclidean distances from current test sample to all training samples
        dist = sum((X_s_training - X_s_test(i, :)).^2, 2);
        % Sort distances and get indices of nearest neighbors
        [~, indices] = sort(dist);
        nn_indices = indices(1:5);  % Get indices of the nearest neighbors

        % Calculate weights for averaging the predictions
        ones_vect=ones(length(X_s_training),1);
        w_i = ones_vect(nn_indices) / 5;

        % Predict the outcome for the test sample
        y_predicted(i) = sum(w_i .* theta_s(nn_indices));
    end

    %% Evaluate Model Performance
    % Compute and store the mean squared error of predictions (using the temporal spatial version of the mse)
    result_rep(1,rep)=mse_new(y_predicted,X_s_test,Weights_test);
end
result_rep

toc  % End timer and display execution time
