

function predicted_categories = nearest_neighbor_classify(train_image_feats, train_labels, test_image_feats)
    model = fitcknn(train_image_feats, train_labels, 'NumNeighbors', 8, 'Distance', 'cosine');
    predicted_categories = predict(model, test_image_feats);
end

