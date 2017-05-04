%% Function to build samples using the entire scene for the logistical regressions analysis
function [ imagery_delta, gis_sample ] = allData( imagery, image_combo, gis_data )
    imagery_size = size(imagery);
    gis_data_size = size(gis_data);
    
    for i = 1:imagery_size(3)
        temp(:,i) = reshape(imagery(:,:,i),(imagery_size(1) * imagery_size(2)),1);
        temp(:,i) = temp(:,i) == 0.25;
    end
    
    x = size(image_combo);
    for i = 1:x(1)
       imagery_delta(:,i) =  temp(:,image_combo(i,1)) ~= 1 & temp(:,image_combo(i,2)) == 1;
    end
    
    for i = 1:gis_data_size(3)
        gis_sample(:,i) = reshape(gis_data(:,:,i),(imagery_size(1) * imagery_size(2)),1);
    end
    %gis_sample(:,i+1) = imagery_delta(:,1) ~= 0.25 & imagery_delta(:,2) == 0.25;
end

