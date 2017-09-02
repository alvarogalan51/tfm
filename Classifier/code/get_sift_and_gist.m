
function sift_and_gists = get_sift_and_gist(image_paths)

    gists = get_gist(image_paths);
    sifts = get_bags_of_sifts(image_paths);
    
    sift_and_gists = [gists sifts];

end