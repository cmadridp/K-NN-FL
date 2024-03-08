function mse = mse_new(y_predicted,X_s_test,Weights_tes)
    diff = y_predicted(:,1) - f_S1_vec_eva(X_s_test);
    mse = sum((diff.^2).*Weights_tes(:,1));
end