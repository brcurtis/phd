%% Setup
clear variables;
clc;

%% Load data
inputPath_gis = '/Users/briancurtis/Documents/PhD/model/data/gis_data/';
inputPath_img = '/Users/briancurtis/Documents/PhD/model/data/imagery/';
images = [ '1990133';'1995163';'2016157' ];
image_combo = [1,2;1,3;2,3];
gis = { 'slope';'expressways';'major_water';'city_centers';'roads';'urban_area';'parks';'schools';'cultural_facilities';'census';'railways';'farmland'; 'industrial' };
[ imagery, gis_data ] = load_data( images, inputPath_img, gis, inputPath_gis );

%% Setup variables and sample set generation

% Determine the size of the imagery
imagery_size = size(imagery);

% Percentage of image size to sample
num_samples = 0.001;

% Generate the samples and prepare the data
%[ imagery_delta, gis_sample ] = allData( imagery, image_combo, gis_data );
[ all_data, sample_set ] = genSamples( imagery, image_combo, gis_data, num_samples );

