function image_feats = get_pyramid_sift(image_paths)

    % vocab = vocab_size x 128
    load('pyramid_vocab.mat')

    % size of histogram or num_features
    d = 400;
    n = length(image_paths);
    
    image_feats = ones(n, d);
    
    % build kdtree for distance querying
    kdtree = vl_kdtreebuild(pyramid_vocab', 'NumTrees', 1);
    
    for i = 1:n
        
        image = single(imread(image_paths{i}));

        descriptors = get_pyramid_descriptors_for_image(image);
        
        [min_indices, ~] = vl_kdtreequery(kdtree, pyramid_vocab', single(descriptors));
        
        % Build histogram from min_indices
        count_histogram = hist(single(min_indices), d);
        
        image_feats(i, :) = count_histogram ./ size(descriptors, 2);
        
        print_progress_string(i, n, 'Getting Pyramid Sifts');
    end

end