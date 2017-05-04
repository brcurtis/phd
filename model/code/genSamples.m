%% Function to generate samples for the logistical regressions analysis
function [ all_data, sample_set ] = genSamples( imagery, image_combo, gis_data, num_samples )    
    
    %% Setup variables
    imagery_size = size(imagery);
    gis_data_size = size(gis_data);
    sample_size = round((imagery_size(1) * imagery_size(2) * num_samples));
    
    %% Build the land use
    for i = 1:imagery_size(3)
        land_use(:,i) = reshape(imagery(:,:,i),(imagery_size(1) * imagery_size(2)),1);
    end
    
    %% Build the urbanization results from the imagery
    for i = 1:imagery_size(3)
        temp(:,i) = reshape(imagery(:,:,i),(imagery_size(1) * imagery_size(2)),1);
        temp(:,i) = temp(:,i) == 22 | temp(:,i) == 23 | temp(:,i) == 24;
    end
    
    x = size(image_combo);
    for i = 1:x(1)
       urb_results(:,i) =  temp(:,image_combo(i,1)) ~= 1 & temp(:,image_combo(i,2)) == 1;
    end
    
    %% Build the GIS library
    for i = 1:gis_data_size(3)
        gis(:,i) = reshape(gis_data(:,:,i),(imagery_size(1) * imagery_size(2)),1);
    end
    
    %% Build the neighborhood datasets
    nhood_size = 3;
    
    % urban area
    omega_urb = nlfilter((imagery(:,:,1) >= 22 & imagery(:,:,1) <= 24),[nhood_size nhood_size],@omega_fun);
    omega_urb = reshape(omega_urb,(imagery_size(1)*imagery_size(2)),1) .* 8;
    
    % non-urban
    % omega_non_urb = 8 - omega_urb;
    % omega_non_urb = nlfilter((imagery(:,:,1) < 22 | imagery(:,:,1) > 24),[nhood_size nhood_size],@omega_fun);
    % omega_non_urb = reshape(omega_non_urb,(imagery_size(1)*imagery_size(2)),1) .* 8;
    
    
    %% Concat all data together
    all_data = cat(2, im2double(land_use), gis, omega_urb, urb_results);
    all_data_size = size(all_data);
    
    %% Generate sample_set
    %land_use_types = unique(land_use);
    %land_use_types_size = size(land_use_types);
    
%     x = size(image_combo);
%     for i = 1:x(1)
    
        temp = all_data(:,all_data_size(2)-2) == 1;
        temp = repmat(temp,1,all_data_size(2));
        temp = temp .* all_data;
        temp(all(temp==0,2),:)=[];
        urb = datasample(temp,round(sample_size/2));

        temp = all_data(:,all_data_size(2)-2) == 0;
        temp = repmat(temp,1,all_data_size(2));
        temp = temp .* all_data;
        temp(all(temp==0,2),:)=[];
        no_urb = datasample(temp,round(sample_size/2));

        sample_set = cat(1,urb,no_urb);
        n = size(sample_set);

        x = randperm(n(1));
        x = x';
        sample_set = sample_set(x,:);
        
%     end

% Convert land use back into original format (multiple by 255)
all_data(:,[1:3]) = all_data(:,[1:3]) * 255;
sample_set(:,[1:3]) = sample_set(:,[1:3]) * 255;


end

%% Code Wasteland

% %% Generate sample_set
%     land_use_types = unique(land_use);
%     land_use_types_size = size(land_use_types);
%     
% %     x = size(image_combo);
% %     for i = 1:x(1)
%     
%         temp = all_data(:,all_data_size(2)-2) == 1;
%         temp = repmat(temp,1,all_data_size(2));
%         temp = temp .* all_data;
%         temp(all(temp==0,2),:)=[];
%         urb = datasample(temp,round(sample_size/2));
% 
%         temp = all_data(:,all_data_size(2)-2) == 0;
%         temp = repmat(temp,1,all_data_size(2));
%         temp = temp .* all_data;
%         temp(all(temp==0,2),:)=[];
%         no_urb = datasample(temp,round(sample_size/2));
% 
%         sample_set = cat(1,urb,no_urb);
%         n = size(sample_set);
% 
%         x = randperm(n(1));
%         x = x';
%         sample_set = sample_set(x,:);
%         
% %     end