%%

%%%%%
%%%%% Region 1 August
%%%%%
fileid=fopen("ys_2_A.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b=cell2mat(a);
%%
% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (-40.5:1:-20.5)';
lon_grid = (50.5:1:110.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b = [];
for aux = 1:n
X_s_b = [X_s_b; matrix_lon_lat];
end

%%
fileid=fopen("AE_Weights_2_A.text");
format="%f";
a=textscan(fileid, format);
Weights_b=cell2mat(a);


%%
n = 170;
m = total;
num_splits = 30;
split_size = floor(m / num_splits);
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits

    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s = X_s_b(index_mask, :);
    y_i_js = y_i_js_b(index_mask, :);
    Weights = Weights_b(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js), round(length(y_i_js)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training = X_s(train_indices, :);
y_i_js_training = y_i_js(train_indices);
Weights_training=Weights(train_indices);
X_s_test = X_s(setdiff(1:length(y_i_js), train_indices), :);
y_i_js_test = y_i_js(setdiff(1:length(y_i_js), train_indices));
Weights_test=Weights(setdiff(1:length(y_i_js), train_indices));
    lambda=[110,100,105,90,95].';
    lambda_opt=exp(15);%cv(X_s_training,Weights_training,y_i_js_training,5,5,lambda,exp(-14.8),5,8,2,2);
    %%
y_predicted = NaN(size(y_i_js_test));
theta_s=admm_knnfl_varying_rho(X_s_training,y_i_js_training,Weights_training,lambda_opt,5,exp(-15),5,8,2,2);
%%
for t = 1:size(y_predicted, 1)
    dist = sum((X_s_training - X_s_test(t, :)).^2, 2);
    dist(ismember(X_s_training, X_s_test(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted(t) = sum(w_i .* theta_s(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted,y_i_js_test,Weights_test);
end




%%
%%%%%
%%%%% Region 1 July
%%%%%
fileid=fopen("ys_2_J.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b_J=cell2mat(a);

% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (-40.5:1:-20.5)';
lon_grid = (50.5:1:110.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b_J = [];
for aux = 1:n
X_s_b_J = [X_s_b_J; matrix_lon_lat];
end

%%
fileid=fopen("AE_Weights_2_J.text");
format="%f";
a=textscan(fileid, format);
Weights_b_J=cell2mat(a);


%%
%%
n = 170;
m = total;
num_splits = 30;
split_size = floor(m / num_splits);
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits
    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b_J, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s_J = X_s_b_J(index_mask, :);
    y_i_js_J = y_i_js_b_J(index_mask, :);
    Weights_J = Weights_b_J(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js_J), round(length(y_i_js_J)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training_J = X_s_J(train_indices, :);
y_i_js_training_J = y_i_js_J(train_indices);
Weights_training_J=Weights_J(train_indices);
X_s_test_J = X_s_J(setdiff(1:length(y_i_js_J), train_indices), :);
y_i_js_test_J = y_i_js_J(setdiff(1:length(y_i_js_J), train_indices));
Weights_test_J=Weights_J(setdiff(1:length(y_i_js_J), train_indices));
    lambda=[exp(-15.2),exp(-20),exp(10.5),exp(15),exp(26)].';
    lambda_opt=exp(-15);%cv(X_s_training,Weights_training,y_i_js_training,5,5,lambda,exp(-14.8),15,8,2,2);
    %%
y_predicted_J = NaN(size(y_i_js_test_J));
theta_s_J=admm_knnfl_varying_rho(X_s_training_J,y_i_js_training_J,Weights_training_J,lambda_opt,5,exp(-15),5,8,2,2);
%%
for t = 1:size(y_predicted_J, 1)
    dist = sum((X_s_training_J - X_s_test_J(t, :)).^2, 2);
    dist(ismember(X_s_training_J, X_s_test_J(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training_J),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted_J(t) = sum(w_i .* theta_s_J(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted_J,y_i_js_test_J,Weights_test_J);
end



%%
%%%%%
%%%%% Region 1 June
%%%%%

fileid=fopen("ys_2_Ju.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b_Ju=cell2mat(a);

% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (-40.5:1:-20.5)';
lon_grid = (50.5:1:110.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b_Ju = [];
for aux = 1:n
X_s_b_Ju = [X_s_b_Ju; matrix_lon_lat];
end

%%
fileid=fopen("AE_Weights_2_Ju.text");
format="%f";
a=textscan(fileid, format);
Weights_b_Ju=cell2mat(a);

%%


%%
n = 170;
m = total;
num_splits = 30;
split_size = floor(m / num_splits);
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits
    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b_Ju, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s_Ju = X_s_b_Ju(index_mask, :);
    y_i_js_Ju = y_i_js_b_Ju(index_mask, :);
    Weights_Ju = Weights_b_Ju(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js_Ju), round(length(y_i_js_Ju)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training_Ju = X_s_Ju(train_indices, :);
y_i_js_training_Ju = y_i_js_Ju(train_indices);
Weights_training_Ju=Weights_Ju(train_indices);
X_s_test_Ju = X_s_Ju(setdiff(1:length(y_i_js_Ju), train_indices), :);
y_i_js_test_Ju = y_i_js_Ju(setdiff(1:length(y_i_js_Ju), train_indices));
Weights_test_Ju=Weights_Ju(setdiff(1:length(y_i_js_Ju), train_indices));
    lambda=[exp(-15.2),exp(-20),exp(10.5),exp(15),exp(26)].';
    lambda_opt=exp(-15);%cv(X_s_training,Weights_training,y_i_js_training,5,5,lambda,exp(-14.8),15,8,2,2);
    %%
y_predicted_Ju = NaN(size(y_i_js_test_Ju));
theta_s_Ju=admm_knnfl_varying_rho(X_s_training_Ju,y_i_js_training_Ju,Weights_training_Ju,lambda_opt,5,exp(-15),5,8,2,2);
%%
for t = 1:size(y_predicted_Ju, 1)
    dist = sum((X_s_training_Ju - X_s_test_Ju(t, :)).^2, 2);
    dist(ismember(X_s_training_Ju, X_s_test_Ju(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training_Ju),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted_Ju(t) = sum(w_i .* theta_s_Ju(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted_Ju,y_i_js_test_Ju,Weights_test_Ju);
end




%%
%%%
%%%august region 2
%%%

fileid=fopen("ys_A_1.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b_A_R2=cell2mat(a);
%%
% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (35.5:1:54.5)';
lon_grid = (165.5:1:224.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b_A_R2 = [];
for aux = 1:n
X_s_b_A_R2 = [X_s_b_A_R2; matrix_lon_lat];
end
%X_s_b = zscore(X_s_b);
%%
fileid=fopen("A_Weights_A_1.text");
format="%f";
a=textscan(fileid, format);
Weights_b_A_R2=cell2mat(a);


%%
n = 170;
m = total;
num_splits = 20;
split_size = floor(m / num_splits);
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits


    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b_A_R2, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s_A_R2 = X_s_b_A_R2(index_mask, :);
    y_i_js_A_R2 = y_i_js_b_A_R2(index_mask, :);
    Weights_A_R2 = Weights_b_A_R2(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js_A_R2), round(length(y_i_js_A_R2)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training_A_R2 = X_s_A_R2(train_indices, :);
y_i_js_training_A_R2 = y_i_js_A_R2(train_indices);
Weights_training_A_R2=Weights_A_R2(train_indices);
X_s_test_A_R2 = X_s_A_R2(setdiff(1:length(y_i_js_A_R2), train_indices), :);
y_i_js_test_A_R2 = y_i_js_A_R2(setdiff(1:length(y_i_js_A_R2), train_indices));
Weights_test_A_R2=Weights_A_R2(setdiff(1:length(y_i_js_A_R2), train_indices));
    lambda=[110,100,105,90,95].';
    lambda_opt=1;%cv(X_s_training,Weights_training,y_i_js_training,5,5,lambda,exp(-14.8),5,8,2,2);
    %%
y_predicted_A_R2 = NaN(size(y_i_js_test_A_R2));
theta_s_A_R2=admm_knnfl_varying_rho(X_s_training_A_R2,y_i_js_training_A_R2,Weights_training_A_R2,lambda_opt,5,exp(-15),5,10,2,2);
%%
for t = 1:size(y_predicted_A_R2, 1)
    dist = sum((X_s_training_A_R2 - X_s_test_A_R2(t, :)).^2, 2);
    dist(ismember(X_s_training_A_R2, X_s_test_A_R2(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training_A_R2),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted_A_R2(t) = sum(w_i .* theta_s_A_R2(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted_A_R2,y_i_js_test_A_R2,Weights_test_A_R2);
end

%%
%%%
%%%july region 2
%%%

fileid=fopen("ys_J_1.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b_J_R2=cell2mat(a);
%%
% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (35.5:1:54.5)';
lon_grid = (165.5:1:224.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b_J_R2 = [];
for aux = 1:n
X_s_b_J_R2 = [X_s_b_J_R2; matrix_lon_lat];
end
%X_s_b = zscore(X_s_b);
%%
fileid=fopen("A_Weights_J_1.text");
format="%f";
a=textscan(fileid, format);
Weights_b_J_R2=cell2mat(a);


%%
n = 170;
m = total;
num_splits = 20;
split_size = floor(m / num_splits);
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits


    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b_J_R2, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s_J_R2 = X_s_b_J_R2(index_mask, :);
    y_i_js_J_R2 = y_i_js_b_J_R2(index_mask, :);
    Weights_J_R2 = Weights_b_J_R2(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js_J_R2), round(length(y_i_js_J_R2)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training_J_R2 = X_s_J_R2(train_indices, :);
y_i_js_training_J_R2 = y_i_js_J_R2(train_indices);
Weights_training_J_R2=Weights_J_R2(train_indices);
X_s_test_J_R2 = X_s_J_R2(setdiff(1:length(y_i_js_J_R2), train_indices), :);
y_i_js_test_J_R2 = y_i_js_J_R2(setdiff(1:length(y_i_js_J_R2), train_indices));
Weights_test_J_R2=Weights_J_R2(setdiff(1:length(y_i_js_J_R2), train_indices));
    lambda=[110,100,105,90,95].';
    lambda_opt=1;%cv(X_s_training,Weights_training,y_i_js_training,5,5,lambda,exp(-14.8),5,8,2,2);
    %%
y_predicted_J_R2 = NaN(size(y_i_js_test_J_R2));
theta_s_J_R2=admm_knnfl_varying_rho(X_s_training_J_R2,y_i_js_training_J_R2,Weights_training_J_R2,lambda_opt,5,exp(-15),5,10,2,2);
%%
for t = 1:size(y_predicted_J_R2, 1)
    dist = sum((X_s_training_J_R2 - X_s_test_J_R2(t, :)).^2, 2);
    dist(ismember(X_s_training_J_R2, X_s_test_J_R2(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training_J_R2),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted_J_R2(t) = sum(w_i .* theta_s_J_R2(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted_J_R2,y_i_js_test_J_R2,Weights_test_J_R2);
end

%%
%%%
%%%june region 2
%%%

fileid=fopen("ys_Ju_1.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b_Ju_R2=cell2mat(a);
%%
% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (35.5:1:54.5)';
lon_grid = (165.5:1:224.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b_Ju_R2 = [];
for aux = 1:n
X_s_b_Ju_R2 = [X_s_b_Ju_R2; matrix_lon_lat];
end
%X_s_b = zscore(X_s_b);
%%
fileid=fopen("A_Weights_Ju_1.text");
format="%f";
a=textscan(fileid, format);
Weights_b_Ju_R2=cell2mat(a);


%%
n = 170;
m = total;
num_splits = 20;
split_size = floor(m / num_splits);
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits


    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b_Ju_R2, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s_Ju_R2 = X_s_b_Ju_R2(index_mask, :);
    y_i_js_Ju_R2 = y_i_js_b_Ju_R2(index_mask, :);
    Weights_Ju_R2 = Weights_b_Ju_R2(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js_Ju_R2), round(length(y_i_js_Ju_R2)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training_Ju_R2 = X_s_Ju_R2(train_indices, :);
y_i_js_training_Ju_R2 = y_i_js_Ju_R2(train_indices);
Weights_training_Ju_R2=Weights_Ju_R2(train_indices);
X_s_test_Ju_R2 = X_s_Ju_R2(setdiff(1:length(y_i_js_Ju_R2), train_indices), :);
y_i_js_test_Ju_R2 = y_i_js_Ju_R2(setdiff(1:length(y_i_js_Ju_R2), train_indices));
Weights_test_Ju_R2=Weights_Ju_R2(setdiff(1:length(y_i_js_Ju_R2), train_indices));
    lambda=[110,100,105,90,95].';
    lambda_opt=exp(10);%cv(X_s_training,Weights_training,y_i_js_training,5,5,lambda,exp(-14.8),5,8,2,2);
    %%
y_predicted_Ju_R2 = NaN(size(y_i_js_test_Ju_R2));
theta_s_Ju_R2=admm_knnfl_varying_rho(X_s_training_Ju_R2,y_i_js_training_Ju_R2,Weights_training_Ju_R2,lambda_opt,5,exp(-15),5,10,2,2);
%%
for t = 1:size(y_predicted_Ju_R2, 1)
    dist = sum((X_s_training_Ju_R2 - X_s_test_Ju_R2(t, :)).^2, 2);
    dist(ismember(X_s_training_Ju_R2, X_s_test_Ju_R2(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training_Ju_R2),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted_Ju_R2(t) = sum(w_i .* theta_s_Ju_R2(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted_Ju_R2,y_i_js_test_Ju_R2,Weights_test_Ju_R2);
end



%%

%%%
%%%august region 3
%%%
fileid=fopen("testA.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b_A_R1=cell2mat(a);
%%
% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (10.5:1:39.5)';
lon_grid = (300.5:1:339.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b_A_R1 = [];
for aux = 1:n
X_s_b_A_R1 = [X_s_b_A_R1; matrix_lon_lat];
end

%%
fileid=fopen("AE_WeightsA.text");
format="%f";
a=textscan(fileid, format);
Weights_b_A_R1=cell2mat(a);


%%
n = 170;
m = total;
num_splits = 10;
split_size = floor(m / num_splits);
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits
    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b_A_R1, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s_A_R1 = X_s_b_A_R1(index_mask, :);
    y_i_js_A_R1 = y_i_js_b_A_R1(index_mask, :);
    Weights_A_R1 = Weights_b_A_R1(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js_A_R1), round(length(y_i_js_A_R1)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training_A_R1 = X_s_A_R1(train_indices, :);
y_i_js_training_A_R1 = y_i_js_A_R1(train_indices);
Weights_training_A_R1=Weights_A_R1(train_indices);
X_s_test_A_R1 = X_s_A_R1(setdiff(1:length(y_i_js_A_R1), train_indices), :);
y_i_js_test_A_R1 = y_i_js_A_R1(setdiff(1:length(y_i_js_A_R1), train_indices));
Weights_test_A_R1=Weights_A_R1(setdiff(1:length(y_i_js_A_R1), train_indices));
    lambda=[exp(-15.2),exp(-10),exp(-10.5),exp(-24),exp(-16)].';
    lambda_opt=exp(15);%cv(X_s_training,ms_training,Weights_training,y_i_js_training,5,5,lambda,0.00001,5,10,2,2);
    %%
y_predicted_A_R1 = NaN(size(y_i_js_test_A_R1));
theta_s_A_R1=admm_knnfl_varying_rho(X_s_training_A_R1,y_i_js_training_A_R1,Weights_training_A_R1,lambda_opt,5,exp(-14.6),5,10,2,2);
%%
for t = 1:size(y_predicted_A_R1, 1)
    dist = sum((X_s_training_A_R1 - X_s_test_A_R1(t, :)).^2, 2);
    dist(ismember(X_s_training_A_R1, X_s_test_A_R1(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training_A_R1),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted_A_R1(t) = sum(w_i .* theta_s_A_R1(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted_A_R1,y_i_js_test_A_R1,Weights_test_A_R1);
end

%%
%%%
%%%july region 3
%%%




fileid=fopen("testJ.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b_J_R1=cell2mat(a);

% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (10.5:1:39.5)';
lon_grid = (300.5:1:339.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b_J_R1 = [];
for aux = 1:n
X_s_b_J_R1 = [X_s_b_J_R1; matrix_lon_lat];
end

%%
fileid=fopen("AE_WeightsJ.text");
format="%f";
a=textscan(fileid, format);
Weights_b_J_R1=cell2mat(a);


%%
n = 170;
m = total;
num_splits = 10;
split_size = floor(m / num_splits)
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits
    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b_J_R1, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s_J_R1 = X_s_b_J_R1(index_mask, :);
    y_i_js_J_R1 = y_i_js_b_J_R1(index_mask, :);
    Weights_J_R1 = Weights_b_J_R1(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js_J_R1), round(length(y_i_js_J_R1)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training_J_R1 = X_s_J_R1(train_indices, :);
y_i_js_training_J_R1 = y_i_js_J_R1(train_indices);
Weights_training_J_R1=Weights_J_R1(train_indices);
X_s_test_J_R1 = X_s_J_R1(setdiff(1:length(y_i_js_J_R1), train_indices), :);
y_i_js_test_J_R1 = y_i_js_J_R1(setdiff(1:length(y_i_js_J_R1), train_indices));
Weights_test_J_R1=Weights_J_R1(setdiff(1:length(y_i_js_J_R1), train_indices));
    lambda=[exp(-15.2),exp(-10),exp(-10.5),exp(-24),exp(-16)].';
    lambda_opt=exp(15.2);%cv(X_s_training,ms_training,Weights_training,y_i_js_training,5,5,lambda,0.00001,5,10,2,2);
    %%
y_predicted_J_R1 = NaN(size(y_i_js_test_J_R1));
theta_s_J_R1=admm_knnfl_varying_rho(X_s_training_J_R1,y_i_js_training_J_R1,Weights_training_J_R1,lambda_opt,5,exp(-15),5,8,2,2);
%%
for t = 1:size(y_predicted_J_R1, 1)
    dist = sum((X_s_training_J_R1 - X_s_test_J_R1(t, :)).^2, 2);
    dist(ismember(X_s_training_J_R1, X_s_test_J_R1(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training_J_R1),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted_J_R1(t) = sum(w_i .* theta_s_J_R1(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted_J_R1,y_i_js_test_J_R1,Weights_test_J_R1);
end


%%
%%%
%%%june region 3
%%%



fileid=fopen("testJu.text");
format="%f";
a=textscan(fileid, format);
y_i_js_b_Ju_R1=cell2mat(a);

% The region of interest is a 10x5 grid, i.e. n=50
lat_grid = (10.5:1:39.5)';
lon_grid = (300.5:1:339.5)';

total = length(lat_grid)*length(lon_grid);
matrix_lon_lat = NaN(total, 2);

for i = 1:length(lat_grid)
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 1) = lon_grid;
matrix_lon_lat(((i-1)*length(lon_grid)+1):(i*length(lon_grid)), 2) = lat_grid(i);
end

n = 170;
% the responses lies in dimension 3
m = total;
% the observations lies in dimension 1
d = 2;
% we define the n_ts
nts = repmat(m, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
% creating X_s
X_s_b_Ju_R1 = [];
for aux = 1:n
X_s_b_Ju_R1 = [X_s_b_Ju_R1; matrix_lon_lat];
end

%%
fileid=fopen("AE_WeightsJu.text");
format="%f";
a=textscan(fileid, format);
Weights_b_Ju_R1=cell2mat(a);


%%
n = 170;
m = total;
num_splits = 10;
split_size = floor(m / num_splits);
nts = repmat(split_size, n, 1);
nm = sum(nts);
disp(nm);
ms=nts;
result_rep=zeros(num_splits,1);
for s = 1:num_splits
    % get the current split indices
    start_index = ((s-1)*split_size) + 1;
    end_index = s * split_size;
    index_mask = false(size(X_s_b_Ju_R1, 1), 1);
    % create index mask for the data within the current split
    for j = 1:n
        index_mask((start_index+(j-1)*m):(end_index+(j-1)*m)) = true(split_size, 1);
    end
    % get the data for the current split
    X_s_Ju_R1 = X_s_b_Ju_R1(index_mask, :);
    y_i_js_Ju_R1 = y_i_js_b_Ju_R1(index_mask, :);
    Weights_Ju_R1 = Weights_b_Ju_R1(index_mask, :);

     %%creating the data for the CV

    % create the training and test indices
train_indices = datasample(1:length(y_i_js_Ju_R1), round(length(y_i_js_Ju_R1)*0.75), 'Replace', false);

% subset the data into training and test sets
X_s_training_Ju_R1 = X_s_Ju_R1(train_indices, :);
y_i_js_training_Ju_R1 = y_i_js_Ju_R1(train_indices);
Weights_training_Ju_R1=Weights_Ju_R1(train_indices);
X_s_test_Ju_R1 = X_s_Ju_R1(setdiff(1:length(y_i_js_Ju_R1), train_indices), :);
y_i_js_test_Ju_R1 = y_i_js_Ju_R1(setdiff(1:length(y_i_js_Ju_R1), train_indices));
Weights_test_Ju_R1=Weights_Ju_R1(setdiff(1:length(y_i_js_Ju_R1), train_indices));
    lambda=[exp(-15.2),exp(-10),exp(-10.5),exp(-24),exp(-16)].';
    lambda_opt=exp(5);%cv(X_s_training,ms_training,Weights_training,y_i_js_training,5,5,lambda,0.00001,5,10,2,2);
    %%
y_predicted_Ju_R1 = NaN(size(y_i_js_test_Ju_R1));
theta_s_Ju_R1=admm_knnfl_varying_rho(X_s_training_Ju_R1,y_i_js_training_Ju_R1,Weights_training_Ju_R1,lambda_opt,5,exp(-14.9),5,10,2,2);
%%
for t = 1:size(y_predicted_Ju_R1, 1)
    dist = sum((X_s_training_Ju_R1 - X_s_test_Ju_R1(t, :)).^2, 2);
    dist(ismember(X_s_training_Ju_R1, X_s_test_Ju_R1(t,:), 'rows')) = inf;

    [~, indices] = sort(dist);
    nn_indices = indices(1:5);
    ones_vect=ones(length(X_s_training_Ju_R1),1);
    w_i = ones_vect(nn_indices) / 5;
    y_predicted_Ju_R1(t) = sum(w_i .* theta_s_Ju_R1(nn_indices));
end
result_rep(s,1)=mse_new(y_predicted_Ju_R1,y_i_js_test_Ju_R1,Weights_test_Ju_R1);
end

%%

%%%%
%%%%Figure KNNFL (all Regions)
%%%%



% Define the range of x and y values
x_range = linspace(165, 225, 100);
y_range = linspace(35, 55, 100);

% Create a grid of x and y values
[X, Y] = meshgrid(x_range, y_range);
x = [X(:), Y(:)];

% Evaluate the function for each point in the grid
Z_A = NaN(size(X));
for i = 1:numel(X)
    y_predicted_A_R2 = y_prediction(x(i,:), X_s_training_A_R2, theta_s_A_R2);
    Z_A(i) = y_predicted_A_R2;
end
Z_A = reshape(Z_A, size(X));

% Evaluate the function for each point in the grid
Z_J = NaN(size(X));
for i = 1:numel(X)
    y_predicted_J_R2 = y_prediction(x(i,:), X_s_training_J_R2, theta_s_J_R2);
    Z_J(i) = y_predicted_J_R2;
end
Z_J = reshape(Z_J, size(X));

% Evaluate the function for each point in the grid
Z_Ju = NaN(size(X));
for i = 1:numel(X)
    y_predicted_Ju_R2 = y_prediction(x(i,:), X_s_training_Ju_R2, theta_s_Ju_R2);
    Z_Ju(i) = y_predicted_Ju_R2;
end
Z_Ju = reshape(Z_Ju, size(X));

% Get the min and max values of Z across all three heatmaps
z_min = min(Z_A(:));
z_max = max(Z_A(:));

% Create a heatmap of the function values for each dataset
figure('Name', 'Heatmaps of Predicted Temperature for Different Months and Locations');
subplot(3,3,6)
imagesc(x_range, y_range, Z_A, [z_min, z_max]);
title('August{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_J = min(Z_J(:));
z_max_J = max(Z_J(:));

subplot(3,3,5)
imagesc(x_range, y_range, Z_J, [z_min_J, z_max_J]);
title('July{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_Ju = min(Z_Ju(:));
z_max_Ju = max(Z_Ju(:));
subplot(3,3,4)
imagesc(x_range, y_range, Z_Ju, [z_min_Ju, z_max_Ju]);
title('June{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


% Define the range of x and y values
x_range = linspace(300, 350, 100);
y_range = linspace(10, 40, 100);
% Create a grid of x and y values
[X, Y] = meshgrid(x_range, y_range);
x = [X(:), Y(:)];

% Evaluate the function for each point in the grid
Z_A = NaN(size(X));
for i = 1:numel(X)
    y_predicted_A_R1 = y_prediction(x(i,:), X_s_training_A_R1, theta_s_A_R1);
    Z_A(i) = y_predicted_A_R1;
end
Z_A = reshape(Z_A, size(X));

% Evaluate the function for each point in the grid
Z_J = NaN(size(X));
for i = 1:numel(X)
    y_predicted_J_R1 = y_prediction(x(i,:), X_s_training_J_R1, theta_s_J_R1);
    Z_J(i) = y_predicted_J_R1;
end
Z_J = reshape(Z_J, size(X));

% Evaluate the function for each point in the grid
Z_Ju = NaN(size(X));
for i = 1:numel(X)
    y_predicted_Ju_R1 = y_prediction(x(i,:), X_s_training_Ju_R1, theta_s_Ju_R1);
    Z_Ju(i) = y_predicted_Ju_R1;
end
Z_Ju = reshape(Z_Ju, size(X));

% Get the min and max values of Z across all three heatmaps
z_min = min(Z_A(:));
z_max = max(Z_A(:));


subplot(3, 3, 3);
imagesc(x_range, y_range, Z_A, [z_min, z_max]);
title('August{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_J = min(Z_J(:));
z_max_J = max(Z_J(:));

subplot(3, 3, 2);
imagesc(x_range, y_range, Z_J, [z_min_J, z_max_J]);
title('July{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_Ju = min(Z_Ju(:));
z_max_Ju = max(Z_Ju(:));
subplot(3, 3, 1);
imagesc(x_range, y_range, Z_Ju, [z_min_Ju, z_max_Ju]);
title('June{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;




% Define the range of x and y values
x_range = linspace(50, 110, 100);
y_range = linspace(-40, -20, 100);

% Create a grid of x and y values
[X, Y] = meshgrid(x_range, y_range);
x = [X(:), Y(:)];

% Evaluate the function for each point in the grid
Z_A = NaN(size(X));
for i = 1:numel(X)
    y_predicted_A = y_prediction(x(i,:), X_s_training, theta_s);
    Z_A(i) = y_predicted_A;
end
Z_A = reshape(Z_A, size(X));

% Evaluate the function for each point in the grid
Z_J = NaN(size(X));
for i = 1:numel(X)
    y_predicted_J = y_prediction(x(i,:), X_s_training_J, theta_s_J);
    Z_J(i) = y_predicted_J;
end
Z_J = reshape(Z_J, size(X));

% Evaluate the function for each point in the grid
Z_Ju = NaN(size(X));
for i = 1:numel(X)
    y_predicted_Ju = y_prediction(x(i,:), X_s_training_Ju, theta_s_Ju);
    Z_Ju(i) = y_predicted_Ju;
end
Z_Ju = reshape(Z_Ju, size(X));

% Get the min and max values of Z across all three heatmaps
z_min = min(Z_A(:));
z_max = max(Z_A(:));

% Create a heatmap of the function values for each dataset
subplot(3,3,9)
imagesc(x_range, y_range, Z_A, [z_min, z_max]);
title('August{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_J = min(Z_J(:));
z_max_J = max(Z_J(:));

subplot(3,3,8)
imagesc(x_range, y_range, Z_J, [z_min_J, z_max_J]);
title('July{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_Ju = min(Z_Ju(:));
z_max_Ju = max(Z_Ju(:));
subplot(3,3,7)
imagesc(x_range, y_range, Z_Ju, [z_min_Ju, z_max_Ju]);
title('June Admm-Knnfl');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


%%
%%%%
%%%%Figure MARS vs KNNFL (Region 2)
%%%%


% Define the range of x and y values
x_range = linspace(165, 225, 100);
y_range = linspace(35, 55, 100);

% Create a grid of x and y values
[X, Y] = meshgrid(x_range, y_range);
x = [X(:), Y(:)];

% Evaluate the function for each point in the grid
Z_A = NaN(size(X));
for i = 1:numel(X)
    y_predicted_A_R2 = y_prediction(x(i,:), X_s_training_A_R2, theta_s_A_R2);
    Z_A(i) = y_predicted_A_R2;
end
Z_A = reshape(Z_A, size(X));

% Evaluate the function for each point in the grid
Z_J = NaN(size(X));
for i = 1:numel(X)
    y_predicted_J_R2 = y_prediction(x(i,:), X_s_training_J_R2, theta_s_J_R2);
    Z_J(i) = y_predicted_J_R2;
end
Z_J = reshape(Z_J, size(X));

% Evaluate the function for each point in the grid
Z_Ju = NaN(size(X));
for i = 1:numel(X)
    y_predicted_Ju_R2 = y_prediction(x(i,:), X_s_training_Ju_R2, theta_s_Ju_R2);
    Z_Ju(i) = y_predicted_Ju_R2;
end
Z_Ju = reshape(Z_Ju, size(X));

% Get the min and max values of Z across all three heatmaps
z_min = min(Z_A(:));
z_max = max(Z_A(:));

% Create a heatmap of the function values for each dataset
figure('Name', 'Heatmaps of Predicted Temperature for Different Months and Locations');
subplot(2,3,3)
imagesc(x_range, y_range, Z_A, [z_min, z_max]);
title('August{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_J = min(Z_J(:));
z_max_J = max(Z_J(:));

subplot(2,3,2)
imagesc(x_range, y_range, Z_J, [z_min_J, z_max_J]);
title('July{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_Ju = min(Z_Ju(:));
z_max_Ju = max(Z_Ju(:));
subplot(2,3,1)
imagesc(x_range, y_range, Z_Ju, [z_min_Ju, z_max_Ju]);
title('June{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


% Read data from files
fileid=fopen("marsA_1.text");
format="%f";
a=textscan(fileid, format);
pred_mars=cell2mat(a);
pred_mars= reshape(pred_mars, size(X));
fclose(fileid);

fileid=fopen("marsJ_1.text");
format="%f";
a=textscan(fileid, format);
pred_mars_J=cell2mat(a);
pred_mars_J= reshape(pred_mars_J, size(X));
fclose(fileid);

fileid=fopen("marsJu_1.text");
format="%f";
a=textscan(fileid, format);
pred_mars_Ju=cell2mat(a);
pred_mars_Ju= reshape(pred_mars_Ju, size(X));
fclose(fileid);

% Get the min and max values of Z across all three heatmaps
pred_min = min([pred_mars(:); pred_mars_J(:); pred_mars_Ju(:)]);
pred_max = max([pred_mars(:); pred_mars_J(:); pred_mars_Ju(:)]);

% Create subplots
subplot(2, 3, 6);
imagesc(x_range, y_range, pred_mars,[z_min, z_max]);
title('August MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


subplot(2, 3, 5);
imagesc(x_range, y_range, pred_mars_J, [z_min_J, z_max_J]);
title('July MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


subplot(2, 3, 4);
imagesc(x_range, y_range, pred_mars_Ju, [z_min_Ju, z_max_Ju]);
title('June MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;
set(gcf,'color','w');



%%


%%%%
%%%%Figure MARS vs KNNFL (Region 3)
%%%%


% Define the range of x and y values
x_range = linspace(300, 350, 100);
y_range = linspace(10, 40, 100);
% Create a grid of x and y values
[X, Y] = meshgrid(x_range, y_range);
x = [X(:), Y(:)];

% Evaluate the function for each point in the grid
Z_A = NaN(size(X));
for i = 1:numel(X)
    y_predicted_A_R1 = y_prediction(x(i,:), X_s_training_A_R1, theta_s_A_R1);
    Z_A(i) = y_predicted_A_R1;
end
Z_A = reshape(Z_A, size(X));

% Evaluate the function for each point in the grid
Z_J = NaN(size(X));
for i = 1:numel(X)
    y_predicted_J_R1 = y_prediction(x(i,:), X_s_training_J_R1, theta_s_J_R1);
    Z_J(i) = y_predicted_J_R1;
end
Z_J = reshape(Z_J, size(X));

% Evaluate the function for each point in the grid
Z_Ju = NaN(size(X));
for i = 1:numel(X)
    y_predicted_Ju_R1 = y_prediction(x(i,:), X_s_training_Ju_R1, theta_s_Ju_R1);
    Z_Ju(i) = y_predicted_Ju_R1;
end
Z_Ju = reshape(Z_Ju, size(X));

% Get the min and max values of Z across all three heatmaps
z_min = min(Z_A(:));
z_max = max(Z_A(:));


subplot(2, 3, 3);
imagesc(x_range, y_range, Z_A, [z_min, z_max]);
title('August{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_J = min(Z_J(:));
z_max_J = max(Z_J(:));

subplot(2, 3, 2);
imagesc(x_range, y_range, Z_J, [z_min_J, z_max_J]);
title('July{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_Ju = min(Z_Ju(:));
z_max_Ju = max(Z_Ju(:));
subplot(2, 3, 1);
imagesc(x_range, y_range, Z_Ju, [z_min_Ju, z_max_Ju]);
title('June{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

% Read data from files
fileid=fopen("marsA.text");
format="%f";
a=textscan(fileid, format);
pred_mars=cell2mat(a);
pred_mars= reshape(pred_mars, size(X));
fclose(fileid);

fileid=fopen("marsJ.text");
format="%f";
a=textscan(fileid, format);
pred_mars_J=cell2mat(a);
pred_mars_J= reshape(pred_mars_J, size(X));
fclose(fileid);

fileid=fopen("marsJu.text");
format="%f";
a=textscan(fileid, format);
pred_mars_Ju=cell2mat(a);
pred_mars_Ju= reshape(pred_mars_Ju, size(X));
fclose(fileid);

% Get the min and max values of Z across all three heatmaps
pred_min = min([pred_mars(:); pred_mars_J(:); pred_mars_Ju(:)]);
pred_max = max([pred_mars(:); pred_mars_J(:); pred_mars_Ju(:)]);

% Create subplots
subplot(2, 3, 6);
imagesc(x_range, y_range, pred_mars,[z_min, z_max]);
title('August MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


subplot(2, 3, 5);
imagesc(x_range, y_range, pred_mars_J, [z_min_J, z_max_J]);
title('July MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


subplot(2, 3, 4);
imagesc(x_range, y_range, pred_mars_Ju, [z_min_Ju, z_max_Ju]);
title('June MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;
set(gcf,'color','w');


%%



%%%%
%%%%Figure MARS vs KNNFL (Region 1)
%%%%
% Define the range of x and y values
x_range = linspace(50, 110, 100);
y_range = linspace(-40, -20, 100);

% Create a grid of x and y values
[X, Y] = meshgrid(x_range, y_range);
x = [X(:), Y(:)];

% Evaluate the function for each point in the grid
Z_A = NaN(size(X));
for i = 1:numel(X)
    y_predicted_A = y_prediction(x(i,:), X_s_training, theta_s);
    Z_A(i) = y_predicted_A;
end
Z_A = reshape(Z_A, size(X));

% Evaluate the function for each point in the grid
Z_J = NaN(size(X));
for i = 1:numel(X)
    y_predicted_J = y_prediction(x(i,:), X_s_training_J, theta_s_J);
    Z_J(i) = y_predicted_J;
end
Z_J = reshape(Z_J, size(X));

% Evaluate the function for each point in the grid
Z_Ju = NaN(size(X));
for i = 1:numel(X)
    y_predicted_Ju = y_prediction(x(i,:), X_s_training_Ju, theta_s_Ju);
    Z_Ju(i) = y_predicted_Ju;
end
Z_Ju = reshape(Z_Ju, size(X));

% Get the min and max values of Z across all three heatmaps
z_min = min(Z_A(:));
z_max = max(Z_A(:));

% Create a heatmap of the function values for each dataset
figure('Name', 'Heatmaps of Predicted Temperature for Different Months and Locations');
subplot(2,3,3)
imagesc(x_range, y_range, Z_A, [z_min, z_max]);
title('August{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_J = min(Z_J(:));
z_max_J = max(Z_J(:));

subplot(2,3,2)
imagesc(x_range, y_range, Z_J, [z_min_J, z_max_J]);
title('July{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;

z_min_Ju = min(Z_Ju(:));
z_max_Ju = max(Z_Ju(:));
subplot(2,3,1)
imagesc(x_range, y_range, Z_Ju, [z_min_Ju, z_max_Ju]);
title('June{\it K}-NN-FL');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


% Read data from files
fileid=fopen("marsA_2.text");
format="%f";
a=textscan(fileid, format);
pred_mars=cell2mat(a);
pred_mars= reshape(pred_mars, size(X));
fclose(fileid);

fileid=fopen("marsJ_2.text");
format="%f";
a=textscan(fileid, format);
pred_mars_J=cell2mat(a);
pred_mars_J= reshape(pred_mars_J, size(X));
fclose(fileid);

fileid=fopen("marsJu_2.text");
format="%f";
a=textscan(fileid, format);
pred_mars_Ju=cell2mat(a);
pred_mars_Ju= reshape(pred_mars_Ju, size(X));
fclose(fileid);

% Get the min and max values of Z across all three heatmaps
pred_min = min([pred_mars(:); pred_mars_J(:); pred_mars_Ju(:)]);
pred_max = max([pred_mars(:); pred_mars_J(:); pred_mars_Ju(:)]);

% Create subplots
subplot(2, 3, 6);
imagesc(x_range, y_range, pred_mars,[z_min, z_max]);
title('August MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


subplot(2, 3, 5);
imagesc(x_range, y_range, pred_mars_J, [z_min_J, z_max_J]);
title('July MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;


subplot(2, 3, 4);
imagesc(x_range, y_range, pred_mars_Ju, [z_min_Ju, z_max_Ju]);
title('June MARS');
xlabel('Longitude');
ylabel('Latitude');
colorbar;
set(gcf,'color','w');


