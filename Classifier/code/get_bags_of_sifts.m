
function image_feats = get_bags_of_sifts(image_paths)
% image_paths is an N x 1 cell array of strings where each string is an
% image path on the file system.

    % vocab = vocab_size x 128
    load('vocab.mat')

    % size of histogram or num_features
    d = 400;
    n = length(image_paths);
    
    image_feats = ones(n, d);
    
    % build kdtree for distance querying
    kdtree = vl_kdtreebuild(vocab', 'NumTrees', 1);
    
    for i = 1:n
        
        image = single(imread(image_paths{i}));
        
        % Get the locations and descriptors for the image
        % with a fast dense sift with fixed size and step
        %
        % descriptors = 128 x num_features_found
        % locations = 2 x num_features_found
        [~, descriptors] = vl_dsift(image, 'step', 8, 'size', 4, 'fast');
        
        % Find the closest vocab given this image's descriptors
        % each descriptor is found closest to the distance in the vocab
        % matrix
        % min finds the closes of the descriptor to the vocab
        % min_indices = 1 x num_features_found with values from
        % 1..vocab_size
        
        rand_descriptors = descriptors;
        
        [min_indices, ~] = vl_kdtreequery(kdtree, vocab', single(rand_descriptors));
        
        % Build histogram from min_indices
        count_histogram = hist(single(min_indices), d);
        
        image_feats(i, :) = count_histogram ./ size(rand_descriptors, 2);
    end
end



