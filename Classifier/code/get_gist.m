
function gists = get_gist(image_paths)

    n = length(image_paths);

    param.imageSize = [256 256]; % it works also with non-square images
    param.orientationsPerScale = [8 8 8 8];
    param.numberBlocks = 4;
    param.fc_prefilt = 4;

    clear gists
    
    for i = 1:n

        img = imread(image_paths{i});
        [gist, param] = LMgist(img, '', param);
        
        if ~exist('gists', 'var')
            gists = zeros(n, size(gist,2));    
        end

        % feature scaling
        mean_gist = mean(gist);
        min_gist = min(gist);
        max_gist = max(gist);
        
        normalized_gist = (gist - mean_gist) ./ (max_gist - min_gist);
        
        gists(i, :) = normalized_gist;

        print_progress_string(i, n, 'Getting Gist');
    end

end