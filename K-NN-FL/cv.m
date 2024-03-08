function [lambda_opt] = cv(X_s_training, ms_training, Weights_training, y_i_js_training, folds_number, number_of_neighb, lambda, rho, max_iter,mu,tau_inc,tau_dec)
    d = size(X_s_training, 2);
    n=size(ms_training);
    ms_train_test = ms_training / folds_number;
    nm_train_test = sum(ms_train_test);

    ms_train_train = 4 * ms_training / folds_number;
    nm_train_train = sum(ms_train_train);

    new_lambda = sort(lambda);
    mse_lambdas = NaN(folds_number, length(new_lambda));
    lambda_min = NaN(1, length(new_lambda));
    for k = 1:folds_number
        X_s_training_training = NaN(nm_train_train, d);
        y_i_js_training_training = NaN(nm_train_train, 1);
        Weights_training_training = NaN(nm_train_train, 1);
    
        X_s_training_test = NaN(nm_train_test, d);
        y_i_js_training_test = NaN(nm_train_test, 1);
        Weights_training_test = NaN(nm_train_test, 1);
    
        aux_2 = 0;
        aux_1 = 0;
        aux_3 = 0;
        aux_4 = 0;
        aux_5 = 0;
        aux_6 = 0;
        for i = 1:n
            aux = ms_training(i) / folds_number * (k - 1);
            aux_0 = ms_training(i) / folds_number * (k - 1);
            if aux > 0
                X_s_training_training((aux_2+1):(aux_2+aux_0), :) = X_s_training((aux_1+1):(aux_1+(k-1)*ms_training(i)/folds_number), :);
                y_i_js_training_training((aux_2+1):(aux_2+aux_0), :) = y_i_js_training((aux_1+1):(aux_1+(k-1)*ms_training(i)/folds_number), :);
                Weights_training_training((aux_2+1):(aux_2+aux_0), :) = Weights_training((aux_1+1):(aux_1+(k-1)*ms_training(i)/folds_number), :);

                X_s_training_test((aux_5+1):(ms_training(i)/(folds_number)+aux_5), :) = X_s_training((aux_4+(k-1)*ms_training(i)/folds_number+1):(aux_4+k*ms_training(i)/folds_number), :);
                y_i_js_training_test((aux_5+1):(ms_training(i)/(folds_number)+aux_5), :) = y_i_js_training((aux_4+(k-1)*ms_training(i)/folds_number+1):(aux_4+k*ms_training(i)/folds_number), :);
                Weights_training_test((aux_5+1):(ms_training(i)/(folds_number)+aux_5), :) = Weights_training((aux_4+(k-1)*ms_training(i)/folds_number+1):(aux_4+k*ms_training(i)/folds_number), :);
                if aux+(ms_training(i)/folds_number)+1 <= ms_training(i)
                    X_s_training_training(aux_2+aux_0+1:ms_train_train(i)+aux_3,:) = X_s_training(aux_4+k*(ms_training(i)/folds_number)+1:ms_training(i)+aux_4,:);
                    y_i_js_training_training(aux_2+aux_0+1:ms_train_train(i)+aux_3,:) = y_i_js_training(aux_4+k*(ms_training(i)/folds_number)+1:ms_training(i)+aux_4,:);
                    Weights_training_training(aux_2+aux_0+1:ms_train_train(i)+aux_3,:) = Weights_training(aux_4+k*(ms_training(i)/folds_number)+1:ms_training(i)+aux_4,:);
                end
                aux_2=aux_2+ms_train_train(i);
                aux_1=aux_1+ms_training(i);
                aux_5=aux_5+ms_train_test(i);
            else
                X_s_training_test((aux_5+aux_0+1):((k)*(ms_training(i))/(folds_number)+aux_5),:) = X_s_training((aux_4+1):(aux_4+(k)*(ms_training(i))/(folds_number)),:);
                y_i_js_training_test((aux_5+aux_0+1):((k)*(ms_training(i))/(folds_number)+aux_5),:) = y_i_js_training((aux_4+1):(aux_4+(k)*(ms_training(i))/(folds_number)),:);
                Weights_training_test((aux_5+aux_0+1):((k)*(ms_training(i))/(folds_number)+aux_5),:) = Weights_training((aux_4+1):(aux_4+(k)*(ms_training(i))/(folds_number)),:);
        
                X_s_training_training((aux_2+aux_0+1):(ms_train_train(i)+aux_3),:) = X_s_training((aux_4+(k)*(ms_training(i))/(folds_number)+1):(ms_training(i)+aux_4),:);
                y_i_js_training_training((aux_2+aux_0+1):(ms_train_train(i)+aux_3),:) = y_i_js_training((aux_4+(k)*(ms_training(i))/(folds_number)+1):(ms_training(i)+aux_4),:);
                Weights_training_training((aux_2+aux_0+1):(ms_train_train(i)+aux_3),:) = Weights_training((aux_4+(k)*(ms_training(i))/(folds_number)+1):(ms_training(i)+aux_4),:);
                aux_2 = aux_2 + ms_train_train(i);
                aux_5 = aux_5 + ms_train_test(i);
                aux_1 = aux_1 + ms_training(i);
            end
            aux_4 = aux_4 + ms_training(i);
            aux_3 = aux_3 + ms_train_train(i);
            aux_6 = aux_6 + ms_train_test(i);
        end
       
        for j = 1:length(lambda)
            lambda_j = lambda(j);
            y_predicted = NaN(size(y_i_js_training_test));
            theta_s=admm_knnfl_varying_rho(X_s_training_training,y_i_js_training_training,Weights_training_training,lambda_j, number_of_neighb,rho, max_iter,mu,tau_inc,tau_dec);
            for i = 1:size(y_predicted, 1)
                
                dist = sum((X_s_training_training - X_s_training_test(i, :)).^2, 2);
                [~, indices] = sort(dist);
                nn_indices = indices(1:number_of_neighb);
                ones_vect=ones(length(X_s_training_training),1);
                w_i = ones_vect(nn_indices) / number_of_neighb;
                y_predicted(i) = sum(w_i .* theta_s(nn_indices));
            end

            mse_lambdas(k, j) = mse_new(y_predicted,X_s_training_test,Weights_training_test);
        end
    end

    mean_mse_lambdas = mean(mse_lambdas, 1);
    lambda_opt_index = find(mean_mse_lambdas == min(mean_mse_lambdas), 1);
    lambda_opt = new_lambda(lambda_opt_index);
 end
