function [ imagery, gis_data, T ] = load_data( images, inputPath_img, gis, inputPath_gis )

    % Load imagery and perform urbanization analysis
    % Default input path for the imagery
    inputPath = inputPath_img;

    % Store imagery into an array for analysis
    % imagery = ['1990133';'1995163';'2016157'];
    % filename_extent = '_classified_aoi.tif';  % used for "orig" data

    % Determine # of images in imagery array
    img_array_size = size(images);

    % Load the scenes (i.e. imagery)
    for i=1:img_array_size(1)
        % Create the file path + name and read in image
        % fname = strcat(inputPath,images(i,:),'/',images(i,:),filename_extent);    % used for "orig" data
        % fname = strcat(inputPath,images(i,:),'/',images(i,:));                    % used for "orig" data
        fname = strcat(inputPath,images(i,:),'.tif');
        imagery(:,:,i) = imread(fname);

        % Normalize the scene data
        %imagery(:,:,i) = (imagery(:,:,i) - min(min(imagery(:,:,i)))) / (max(max(imagery(:,:,i))) - min(min(imagery(:,:,i))));
    end
    
    % Load GIS data
    % Create scene_data 3d-matrix and setup variable
    gis_size = size(gis);
    inputPath = inputPath_gis;
    slope_max = 10;
    
    for i = 1:gis_size(1)
        % Load and prep all GIS support data
        %if strcmp( gis(i), 'slope' ) == 1
            %gis_data(:,:,i) = slope_analysis(inputPath,strcat(char(gis(i)),'.tif'),slope_max);
        if strcmp( gis(i), 'population' ) == 1
            gis_data(:,:,i) = imread(strcat(inputPath,char(gis(i)),'.tif'));
        elseif strcmp( gis(i), 'slope' ) == 1
            clear temp;
            temp = imread(strcat(inputPath,char(gis(i)),'.tif'));
            gis_data(:,:,i) = temp(:,:,1);
        elseif strcmp( gis(i), 'major_water' ) == 1
            gis_data(:,:,i) = bwdist(imagery(:,:,1) == 11);
        else
            gis_data(:,:,i) = bwdist(imread(strcat(inputPath,char(gis(i)),'.tif')));
        end
        % Normalize the data
        gis_data(:,:,i) = (gis_data(:,:,i) - min(min(gis_data(:,:,i)))) / (max(max(gis_data(:,:,i))) - min(min(gis_data(:,:,i))));
    end
    
    % Trim the gis_data to the extend of the imagery (i.e. if there are 0's
    % in the imagery, there doesn't need to be gis_information
    clear temp;
    temp = imagery(:,:,1) ~= 0;
    gis_data(:,:,1) = gis_data(:,:,1) .* temp;
    
end

%% Code Wasteland
