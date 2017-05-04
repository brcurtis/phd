function [ urban, urban_growth ] = urbanization( imagery, img_array_size )
    % Loop through imagery array and calculate the # of urban cells and the
    % growth from the previous image
    img_array_size = size(imagery);
    for i=1:img_array_size(3)
        % Determine which cells are set to 'urban' (1) based on classification
        scene_urban = imagery(:,:,i) == 0.25;

        % Count the # of urban cells
        urban(i) = sum(scene_urban(:));

        % Determine urban growth based on current and previous cell count
        if i == 1
            urban_growth(i) = 0;
        else
            urban_growth(i) = urban(i) / urban(i-1);
        end
    end
end