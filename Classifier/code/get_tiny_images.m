
function image_feats = get_tiny_images(image_paths)
% image_paths is an N x 1 cell array of strings where each string is an

n = length(image_paths);
image_feats = zeros(n, 256);

for i = 1:n
    image = imread(image_paths{i});
    
    % Get the smaller of the size dimensions
    % and create scale on that
    image_size = size(image);
    [~, index] = min(image_size);
    scale = 16 / image_size(index);
    new_image = imresize(image, scale);
    new_image = new_image(1:16, 1:16);
    new_image_features = reshape(new_image, 1, []);
    image_feats(i, :) = new_image_features / 255.0;
   
end
