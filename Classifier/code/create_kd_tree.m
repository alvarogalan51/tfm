
% deprecated use vl_kdtree

% points = d+1 x N
% points is a matrix of points with the first row being the unique index
% for creating tree lookup matrix
%
% e.g [ address_1 x_1 y_1 ; address_2 x_2 y_2 ; ... ; address_n x_n, y_n ]
%
% d = dimension of points
% N number of points
% depth of tree
function tree = create_kd_tree(points, depth)

    d = size(points, 1)-1;
    n = size(points, 2);
        
    % Basecase return single node
    if n == 1
        tree = [points(1, 1) -1];
    elseif n == 0
        tree = [];
    elseif n == 2
        point_1 = points(1, 1);
        point_2 = points(1, 2);
        tree = [point_1 point_2; point_2 -1];
    else
        % Alternate sort on each dimension
        sort_on = mod(depth, d) + 1;
        
        % Sort points and grab middle point
        % Based on sort dimension decided earlier
        sorted_points = sortrows(points', sort_on)';
        middle_index = size(sorted_points, 2)/2;
        middle_index = ceil(middle_index);
        middle_point = sorted_points(:, middle_index);
        middle_address = middle_point(1);
        
        % Split points into left and right
        left_points = sorted_points(:, 1:middle_index-1);
        right_points = sorted_points(:, middle_index+1:end);
        
        % Create left and right trees
        left_tree = create_kd_tree(left_points, depth+1);
        right_tree = create_kd_tree(right_points, depth+1);
        
        % Grab root addresses of left and right trees
        left_address = left_tree(1, 1);
        right_address = right_tree(1, 1);
        
        % Add left and right trees to current tree with 
        % middle address and left and right address appended
        tree = [middle_address left_address ; ...
                     middle_address right_address ; ...
                     left_tree ; ...
                     right_tree
                    ];
        
    end
end