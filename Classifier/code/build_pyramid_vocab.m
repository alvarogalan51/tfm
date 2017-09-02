

function vocab = build_pyramid_vocab( image_paths )
% The inputs are images, a N x 1 cell array of image paths and the size of 
% the vocabulary.
    
    % Sample random permutation of images
    sample_size = length(image_paths);
    perm = randperm(length(image_paths), sample_size);
    
    for i = 1:sample_size
        image_index = perm(i);
        image = single(imread(image_paths{image_index}));
        
        descriptors = get_pyramid_descriptors_for_image(image);
        
        num_descriptors = size(descriptors, 2);
        if ~exist('descriptors', 'var')
            total_descriptors = zeros(size(descriptors, 1), num_descriptors * sample_size);
        end
        
        start_index = (i - 1)*num_descriptors + 1;
        end_index = i*num_descriptors;
        total_descriptors(:, start_index:end_index) = descriptors;
        
        print_progress_string(i, sample_size, 'Creating Pyramid Vocab');
    end

    num_clusters = 200;
    
    fprintf('Creating KMeans clusters (%i clusters)', num_clusters);
    % Cluster results
    % 
    % centers = 128 x num_clusters
    % assignments = 1 x num_features_found
    [centers, ~] = vl_kmeans(single(total_descriptors), num_clusters);
    
    % vocab = vocab_size x 128
    vocab = centers';
end





