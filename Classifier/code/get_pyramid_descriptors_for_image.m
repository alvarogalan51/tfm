function descriptors = get_pyramid_descriptors_for_image(image)

    L = 2;
    l0_d = 256;
    l1_d = 128;
    l2_d = 64;

    % Resize to 256x256
    image = imresizecrop(image, [l0_d, l0_d]);

    % L = [0, 1, 2]

    % Descriptors at lev
    % descriptors = 128 x num_features_found
    % locations = 2 x num_features_found

    % l = 0
    l0_descriptors = get_pyramid_descriptor_for_level(image, 0, l0_d);

    % l = 1
    l1_descriptors = get_pyramid_descriptor_for_level(image, 1, l1_d);
    % subtract duplicate descriptors
    l1_descriptors = setdiff(l1_descriptors', l0_descriptors', 'rows')';
    
    % l = 2
    l2_descriptors = get_pyramid_descriptor_for_level(image, 2, l2_d);
    l2_descriptors = setdiff(l2_descriptors', l1_descriptors', 'rows')';

    % Weight the levels and combine
    descriptors = [l0_descriptors .* 1/4, l1_descriptors .* 1/4, l2_descriptors .* 1/2];
    
end