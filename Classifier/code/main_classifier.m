
% FEATURE = 'tiny image';
% FEATURE = 'bag of sift';
 FEATURE = 'gist';
% FEATURE = 'gist & sift';
% FEATURE = 'fisher';
% FEATURE = 'gist & fisher';
% FEATURE = 'spatial pyramid sift';
% FEATURE = 'saved';
% FEATURE = 'spatial pyramid sift + fisher + gist';

% CLASSIFIER = 'nearest neighbor';
CLASSIFIER = 'support vector machine';

% set up paths to VLFeat functions. 
run('D:\Dropbox\Máster en Ciberseguridad\TFM\Definitivo\Classifier\vlfeat-0.9.20\toolbox\vl_setup')
%run('D:\Dropbox\Dropbox\Máster en Ciberseguridad\TFM\vlfeat-0.9.20\toolbox\vl_setup')

data_path = 'D:\Dropbox\Máster en Ciberseguridad\TFM\Definitivo\Classifier\input';
%data_path = 'D:\Dropbox\Dropbox\Máster en Ciberseguridad\TFM\Classifier';

%categories = {'Kitchen', 'Store', 'Bedroom', 'LivingRoom', 'Office', ...
%       'Industrial', 'Suburb', 'InsideCity', 'TallBuilding', 'Street', ...
%       'Highway', 'OpenCountry', 'Coast', 'Mountain', 'Forest'};
categories = {'m5', 'motog5'};


%This list of shortened category names is used later for visualization.
%abbr_categories = {'Kit', 'Sto', 'Bed', 'Liv', 'Off', 'Ind', 'Sub', ...
%    'Cty', 'Bld', 'St', 'HW', 'OC', 'Cst', 'Mnt', 'For'};
abbr_categories = {'m5', 'motog5'};

num_train_per_cat = 60;
num_test_per_cat = 40;

fprintf('Getting paths and labels for all train and test data\n')
[train_image_paths, test_image_paths, train_labels, test_labels] = ...
    get_image_paths(data_path, categories, num_train_per_cat, num_test_per_cat);
%   train_image_paths  1500x1   cell      
%   test_image_paths   1500x1   cell           
%   train_labels       1500x1   cell         
%   test_labels        1500x1   cell          

fprintf('Using %s representation for images\n', FEATURE)

switch lower(FEATURE) 
    
    case 'spatial pyramid sift'
        
          if ~exist('pyramid_vocab.mat', 'file')
            fprintf('No existing pyramid vocabulary found. Computing one from training images\n')
            pyramid_vocab = build_pyramid_vocab(train_image_paths);
            save('pyramid_vocab.mat', 'pyramid_vocab');
          end
        
        train_image_feats = get_pyramid_sift(train_image_paths);
        test_image_feats  = get_pyramid_sift(test_image_paths);
        
    case 'spatial pyramid sift + fisher + gist'
        
        if ~exist('pyramid_vocab.mat', 'file')
           fprintf('No existing pyramid vocabulary found. Computing one from training images\n')
           pyramid_vocab = build_pyramid_vocab(train_image_paths);
            save('pyramid_vocab.mat', 'pyramid_vocab');
        end

        if ~exist('means.mat', 'file')
            fprintf('No existing fisher word vocabulary found. Computing one from training images\n')
            [means, covariances, priors] = build_fisher_vocabulary(train_image_paths);
            save('means.mat', 'means');
            save('covariances.mat', 'covariances');
            save('priors.mat', 'priors');
        end
        
        train_image_feats = get_pyramid_gist_fisher(train_image_paths);
        test_image_feats  = get_pyramid_gist_fisher(test_image_paths);
        
    case 'tiny image'
        train_image_feats = get_tiny_images(train_image_paths);
        test_image_feats  = get_tiny_images(test_image_paths);
        
    case 'placeholder'
        train_image_feats = [];
        test_image_feats = [];
       
    case 'bag of sift'
        if ~exist('vocab.mat', 'file')
            fprintf('No existing visual word vocabulary found. Computing one from training images\n')
            vocab_size = 800; %Larger values will work better (to a point) but be slower to compute
            vocab = build_vocabulary(train_image_paths, vocab_size);
            save('vocab.mat', 'vocab');
        end
        
        % YOU CODE get_bags_of_sifts.m
        train_image_feats = get_bags_of_sifts(train_image_paths);
        test_image_feats  = get_bags_of_sifts(test_image_paths);
        
    case 'gist'
        train_image_feats = get_gist(train_image_paths);
        test_image_feats = get_gist(test_image_paths);
        
    case 'gist & sift'
        train_image_feats = get_sift_and_gist(train_image_paths);
        test_image_feats = get_sift_and_gist(test_image_paths);
     
    case 'gist & fisher'
        
        if ~exist('means.mat', 'file')
            fprintf('No existing fisher word vocabulary found. Computing one from training images\n')
            [means, covariances, priors] = build_fisher_vocabulary(train_image_paths);
            save('means.mat', 'means');
            save('covariances.mat', 'covariances');
            save('priors.mat', 'priors');
        end

        train_image_feats = get_fisher_and_gist(train_image_paths);
        test_image_feats = get_fisher_and_gist(test_image_paths);
        
    case 'fisher'
        
        if ~exist('means.mat', 'file')
            fprintf('No existing fisher word vocabulary found. Computing one from training images\n')
            [means, covariances, priors] = build_fisher_vocabulary(train_image_paths);
            save('means.mat', 'means');
            save('covariances.mat', 'covariances');
            save('priors.mat', 'priors');
        end
        
        train_image_feats = get_fisher_sifts(train_image_paths);
        test_image_feats = get_fisher_sifts(test_image_paths);
        
    case 'saved'
        
        load('test_image_feats.mat')
        load('train_image_feats.mat')
        
    otherwise
        error('Unknown feature type')
end

fprintf('Using %s classifier to predict test set categories\n', CLASSIFIER)

switch lower(CLASSIFIER)    
    case 'nearest neighbor'
        % YOU CODE nearest_neighbor_classify.m 
        predicted_categories = nearest_neighbor_classify(train_image_feats, train_labels, test_image_feats);
        
    case 'support vector machine'
        % YOU CODE svm_classify.m 
        predicted_categories = svm_classify(train_image_feats, train_labels, test_image_feats);
        
    case 'placeholder'
        %The placeholder classifier simply predicts a random category for
        %every test case
        random_permutation = randperm(length(test_labels));
        predicted_categories = test_labels(random_permutation); 
        
    otherwise
        error('Unknown classifier type')
end

create_results_webpage( train_image_paths, ...
                        test_image_paths, ...
                        train_labels, ...
                        test_labels, ...
                        categories, ...
                        abbr_categories, ...
                        predicted_categories)
