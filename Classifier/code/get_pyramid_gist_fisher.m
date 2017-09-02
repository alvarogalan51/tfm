function image_features = get_pyramid_gist_fisher(image_paths)

    gists = get_gist(image_paths);
    sifts = get_fisher_sifts(image_paths);
    pyramids = get_pyramid_sift(image_paths);
    
    image_features = [gists sifts pyramids];

end