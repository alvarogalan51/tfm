
function image_feats = get_fisher_and_gist(image_paths)

    gists = get_gist(image_paths);
    sifts = get_fisher_sifts(image_paths);
    
    image_feats = [gists sifts];

end