
function vocab = build_vocabulary( image_paths, vocab_size )
% The inputs are images, a N x 1 cell array of image paths and the size of 
% the vocabulary.

    total_descriptors = [];
    
    % Sample random permutation of images
    sample_size = length(image_paths);
    disp(sample_size);
    perm = randperm(length(image_paths), sample_size);
    
    for i = 1:length(perm)
        image_index = perm(i);
        
        image = single(imread(image_paths{image_index}));
        % Get the locations and descriptors for the image
        % with a fast dense sift with fixed size and step
        %
        % descriptors = 128 x num_features_found
        % locations = 2 x num_features_found
        [~, descriptors] = vl_dsift(image, 'step', 8, 'size', 4, 'fast');
        
        total_descriptors = [total_descriptors descriptors];
    end

    num_clusters = vocab_size;
    % Cluster results
    % 
    % centers = 128 x num_clusters
    % assignments = 1 x num_features_found
    [centers, ~] = vl_kmeans(single(total_descriptors), num_clusters);
    
    % vocab = vocab_size x 128
    vocab = centers';
end
