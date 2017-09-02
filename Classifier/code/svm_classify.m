
function predicted_categories = svm_classify(train_image_feats, train_labels, test_image_feats)

unique_labels = unique(train_labels);
num_labels = size(unique_labels, 1);

n = size(train_labels, 1);
m = size(test_image_feats, 1);
predicted_categories_mat = zeros(m, num_labels);

% Create model
for i = 1:num_labels
    label = unique_labels(i);
    one_vs_many = strcmp(label, train_labels) * 2 - 1;
    
    fprintf('Predicting Label %i/%i: %s\n', i, num_labels, char(label));
    
    [W, B] = vl_svmtrain(train_image_feats', one_vs_many, 0.00001);
    predicted = W' * test_image_feats' + B';
    predicted_categories_mat(:, i) = predicted';
    
end

[~, indices] = max(predicted_categories_mat, [], 2);
predicted_categories = unique_labels(indices);

