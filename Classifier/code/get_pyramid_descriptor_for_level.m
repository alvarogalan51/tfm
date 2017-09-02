function descriptors = get_pyramid_descriptor_for_level(image, level, step)

    max_sift_step = 8;
    max_sift_size = 8;
    sift_step = max_sift_step * 2 * 1/pow2(level+1); % 8, 4, 2;
    sift_size = max_sift_size * 2 * 1/pow2(level+1); % 8, 4, 2;

    for i = 0:level
        for j = 0:level

            i_start = step*i + 1;
            i_end = step*(i + 1);
            j_start = step*j + 1;
            j_end = step*(j + 1);


            [~, current_descriptors] = vl_dsift(image(i_start:i_end, j_start:j_end), 'step', sift_step, 'size', sift_size, 'fast');

            num_features = size(current_descriptors, 1);
            num_descriptors = size(current_descriptors, 2);

            if ~exist('descriptors', 'var')
                descriptors = zeros(num_features, num_descriptors*pow2(level));
            end

            % place current increment of descriptors in descriptors
            % matrix
            % 128xd
            current_descriptor_index_start = i*num_descriptors + j*num_descriptors + 1;
            current_descriptor_index_end = current_descriptor_index_start + num_descriptors -1;
            descriptors(:, current_descriptor_index_start:current_descriptor_index_end) = current_descriptors;
        end
    end
end