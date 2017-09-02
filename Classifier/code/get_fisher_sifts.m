% http://cs.brown.edu/courses/cs143/results/proj3/mswang/
% http://www.vlfeat.org/overview/encodings.html

function image_feats = get_fisher_sifts(image_paths)

    load('means.mat');
    load('covariances.mat');
    load('priors.mat');
    
    n = length(image_paths);
    
    clear image_feats;
    for i = 1:n
       
        img = imread(image_paths{i});
        [~, descriptors] = vl_dsift(single(img), 'step', 8, 'size', 4, 'fast');
        
        % dx1
        encoding = vl_fisher(single(descriptors), means, covariances, priors);
        
        if ~exist('image_feats', 'var')
            image_feats = zeros(n, size(encoding, 1));
        end
        
        image_feats(i, :) = encoding';
        
       print_progress_string(i, n, 'Encoding Fisher Sifts');
    end
end


