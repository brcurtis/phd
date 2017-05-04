%% Setup
clear variables;
clc;

tic

%% Load data
inputPath_gis = '/Users/briancurtis/Documents/PhD/model/data/gis_data/';
inputPath_img = '/Users/briancurtis/Documents/PhD/model/data/imagery/';
images = [ '1990133';'1995163';'2016157' ];
image_combo = [1,2;1,3;2,3];
gis = { 'slope';'expressways';'major_water';'city_centers';'roads';'urban_area';'parks';'schools';'cultural_facilities';'census';'railways';'farmland'; 'industrial' };
%gis = { 'major_water';'parks' };
[ imagery, gis_data ] = load_data( images, inputPath_img, gis, inputPath_gis );


%% Urban Growth
[ urban, urban_growth ] = urbanization( imagery );   % Calculate the # of urban cells and the growth from the previous image

%Plot the urban cell count
figure
bar(urban)
set(gca,'xticklabel',imagery)

%% AIS model

Ag_num_samples = 0.0015;    % percentage of image size to sample as a seed for initial library
%Ab_init_size = 0.001;      % perfentage of image size to start the initial anitbody library
Ab_init_size = 400;         % initial anitbody library
Ngen = 150;                 % Defines the number of generations (Ngen) for the AIS model

%[ Ab, Abm, Afm_eval, eval ] = ais_v4( scene_size, Ag_num_samples, Ab_init_size, Ngen, scenes, slope, distancefrom_I25_I40, distancefrom_rio_grande, distancefrom_city_center, distancefrom_roads );
[ Abm, Afm_eval, eval ] = ais_v4( Ag_num_samples, Ab_init_size, Ngen, imagery, image_combo, gis_data );


%% CA model
% for i = 1:1%scene_size(3)
%     [ scene_sim(:,:,i) ] = ca( scenes(:,:,i), Abm(:,[1:6]), slope, distancefrom_I25_I40, distancefrom_rio_grande, distancefrom_city_center );
% end

%% Plot data
% figure
% subplot(2,3,[1 4])
% imshow(scenes(:,:,1),[]);
% 
% subplot(2,3,2)
% imshow(slope,[])
% 
% subplot(2,3,3)
% imshow(distancefrom_city_center,[]);
% 
% subplot(2,3,5)
% imshow(distancefrom_I25_I40,[]);
% 
% subplot(2,3,6)
% imshow(distancefrom_rio_grande,[]);

%% Staging area for obsolete code

