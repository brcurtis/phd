function [slope,avg_slope] = slope_analysis(inputPath,fname,slope_max)

    % Default input path and filename for the srtm data
    %inputPath = '/Users/briancurtis/Documents/PhD/model/data/gis_data/';
    %fname = 'slope_15_data.tif';
    

    % Load slope data
    slope = imread(strcat(inputPath,fname));
    
    % Filter slope data
    %slope_check = slope < 5;
    %slope = slope .* slope_check;
    
    % Load 2016 imagery for slope analysis vs. urban cells
    inputPath = '/Users/briancurtis/Documents/PhD/model/data/imagery/';
    imagery = '2016157';
    filename_extent = '_classified_aoi.tif';
    fname = strcat(inputPath,imagery,'/',imagery,filename_extent);
    scene = imread(fname,'TIF');
    
    % Set scene only to be equal to those that are classified as urban
    scene = scene == 1;     % in classification process, urban = 1
    
    % Eliminate those values in the slope that are urban
    filtered_slope = slope .* scene;
    slope_tot = filtered_slope ~= 0;
    slope_delta = filtered_slope>0 & filtered_slope<slope_max;
    
    avg_slope = sum(slope_delta(:)) / sum(slope_tot(:));
    
end