function A_new = average_assemblies(A)

% Reads the map over assemblies and make make all nodes belonging to each assembly equal to the average over the radial cross section of the assembly for each axial slice
nz = size(A,3);
% Process each z-layer

if nz == 26
    load("Assembly_map.mat","Assembly_map")
    % Create a copy of A to store the results
    A_new = A;
    step= (nz-2)/2;
    for z = 2:step:nz-1
        % Get the current z-layer
        current_layer = A(:,:,z:z+step-1);
        current_layer_map = Assembly_map(:,:,z:z+step-1);

        % Find all unique region numbers in this layer
        region_ids = unique(current_layer_map);

        % Process each region
        for region = region_ids'
            % Find all indices for this region in this layer
            region_mask = (current_layer_map == region);
            % Calculate the average of A values in this region
            region_values = current_layer(region_mask);
            region_avg = mean(region_values);
            % Set all values in this region to the average
            A_new(:,:,z:z+step-1) = A_new(:,:,z:z+step-1) .* (~region_mask) + region_avg * region_mask;
        end
    end

    for z = [1,nz]
        current_layer = A(:,:,z);
        current_layer_map = Assembly_map(:,:,z);
        % Find all indices for this region in this layer
        region = -1;
        region_mask = (current_layer_map == region);
        % Calculate the average of A values in this region
        region_values = current_layer(region_mask);
        region_avg = mean(region_values);
        % Set all values in this region to the average
        A_new(:,:,z) = A_new(:,:,z) .* (~region_mask) + region_avg * region_mask;
    end
else
    load("Assembly_map.mat","Assembly_map_refined")
    Assembly_map = Assembly_map_refined;
    % Create a copy of A to store the results
    A_new = A;
    step = (nz-4)/4; % split the core height into four parts
    for z = 3:step:nz-2
        % Get the current z-layer
        current_layer = A(:,:,z:z+step-1);
        current_layer_map = Assembly_map(:,:,z:z+step-1);

        % Find all unique region numbers in this layer
        region_ids = unique(current_layer_map);

        % Process each region
        for region = region_ids'
            % Find all indices for this region in this layer
            region_mask = (current_layer_map == region);
            % Calculate the average of A values in this region
            region_values = current_layer(region_mask);
            region_avg = mean(region_values);
            % Set all values in this region to the average
            A_new(:,:,z:z+step-1) = A_new(:,:,z:z+step-1) .* (~region_mask) + region_avg * region_mask;
        end
    end
    for z = [1,nz-1]
        current_layer = A(:,:,z:z+1);
        current_layer_map = Assembly_map(:,:,z:z+1);
        % Find all indices for this region in this layer
        region = -1;
        region_mask = (current_layer_map == region);
        % Calculate the average of A values in this region
        region_values = current_layer(region_mask);
        region_avg = mean(region_values);
        % Set all values in this region to the average
        A_new(:,:,z:z+1) = A_new(:,:,z:z+1) .* (~region_mask) + region_avg * region_mask;
    end
end


