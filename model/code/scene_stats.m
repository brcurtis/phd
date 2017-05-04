%% Setup
clear variables;
clc;

%% Load data
inputPath_gis = '/Users/briancurtis/Documents/PhD/model/data/gis_data/';
inputPath_img = '/Users/briancurtis/Documents/PhD/model/data/imagery/';
images = [ '1990133';'1995163';'2016157' ];
image_combo = [1,2;1,3;2,3];
gis = { 'slope';'expressways';'major_water';'city_centers';'roads';'urban_area';'parks';'schools';'cultural_facilities';'census';'railways';'farmland'; 'industrial' };
%gis = { 'major_water';'parks' };
[ imagery, gis_data ] = load_data( images, inputPath_img, gis, inputPath_gis );

%% Setup variables and sample set generation
imagery_size = size(imagery);
[ imagery_delta, gis_sample ] = allData( imagery, image_combo, gis_data );
temp = reshape(imagery(:,:,1),(imagery_size(1) * imagery_size(2)),1);    % add landuse
gis_sample = cat(2,gis_sample,temp);
gis_sample_size = size(gis_sample);

for i = 1:imagery_size(3)
    %num_samples = 0.002;
    %sample_set = genSamples( imagery(:,:,[1,2]), gis_data, num_samples );   % a generated sample set

%% Logistical analysis
    % perform the logistical regression
    [B,dev,stats] = glmfit(gis_sample(:,[1:gis_sample_size(2)]),imagery_delta(:,i),'binomial','link','logit');
    mdl = fitglm(gis_sample(:,[1:gis_sample_size(2)]),imagery_delta(:,i),'Distribution','binomial','link','logit');
    scores = mdl.Fitted.Probability;
    [X,Y,T,AUC(i)] = perfcurve(imagery_delta(:,i),scores,1);
    beta(i,:) = stats.beta;
    
    % plot the data
    figure;
    subplot(1,2,1)
    plot(X,Y)
    hold on
    xlabel('False positive rate')
    ylabel('True positive rate')
    title('ROC for Classification by Logistic Regression')
    fprintf('Iteration %d of %d completed\n',i,imagery_size(3));

    % Build a probability image
    B = B';
    B = repmat(B,gis_sample_size(1),1);
    prob_image = 1 ./ (1 + exp(-(B(:,1) + sum(gis_sample .* B(:,2:end),2))));
    prob_image = reshape(prob_image,imagery_size(1),imagery_size(2));

    % plot the probability map
    colormap('jet');
    subplot(1,2,2); imagesc(prob_image);
    colorbar;
end

% Plot data
%all_urban( scene_data, scene_data_size );
%soil_sand(scene_data, scene_data_size );

%% Scene Analysis
function [ all_data, resp ] = scene_analysis( scene_data, scene_data_size );
    urban = ((scene_data(:,1) ~= 0.25 | scene_data(:,1) == 0.25) & scene_data(:,2) == 0.25);
    non_urban = (scene_data(:,1) ~= 0.25 & scene_data(:,2) ~= 0.25);

    urban = repmat(urban,1,scene_data_size(2));
    non_urban = repmat(non_urban,1,scene_data_size(2));

    urban = urban .* scene_data;
    urban(all(urban==0,2),:) = [];
    urban(:,end+1) = 1;
    non_urban = non_urban .* scene_data;
    non_urban(all(non_urban==0,2),:) = [];
    non_urban(:,end+1) = 0;
    
    all_data = [urban;non_urban];
    resp = all_data(:,9);
    all_data(:,9) = [];
    
end

%% All Urbanization
function [] = all_urban ( scene_data, scene_data_size );
    urban = ((scene_data(:,1) ~= 0.25 | scene_data(:,1) == 0.25) & scene_data(:,2) == 0.25);
    non_urban = (scene_data(:,1) ~= 0.25 & scene_data(:,2) ~= 0.25);

    urban = repmat(urban,1,scene_data_size(2));
    non_urban = repmat(non_urban,1,scene_data_size(2));

    urban = urban .* scene_data;
    urban(all(urban==0,2),:) = [];
    non_urban = non_urban .* scene_data;
    non_urban(all(non_urban==0,2),:) = [];

    num_bins = 10;

    figure('name','Urbanization of all data');
    subplot(1,7,1); histogram(scene_data(:,1),4); hold on; histogram(scene_data(:,2),4); hold off; xlabel('Land Useage')                        % Land Use
    legend('urban','non-urban','Location','northwest');
    subplot(1,7,2); histogram(urban(:,3),num_bins); hold on; histogram(non_urban(:,3),num_bins); hold off; xlabel('Slope')                      % Slope
    subplot(1,7,3); histogram(urban(:,4),num_bins); hold on; histogram(non_urban(:,4),num_bins); hold off; xlabel('Distance from I-25 & I-40')	% Distance from I25/I40
    subplot(1,7,4); histogram(urban(:,5),num_bins); hold on; histogram(non_urban(:,5),num_bins); hold off; xlabel('Distance from Rio Grande')	% Distance from Rio Grande
    subplot(1,7,5); histogram(urban(:,6),num_bins); hold on; histogram(non_urban(:,6),num_bins); hold off; xlabel('Distance from City Center')	% Distance from City Center
    subplot(1,7,6); histogram(urban(:,7),num_bins); hold on; histogram(non_urban(:,7),num_bins); hold off; xlabel('Distance from Roads')        % Distance from Roads
    subplot(1,7,7); histogram(urban(:,8),num_bins); hold on; histogram(non_urban(:,8),num_bins); hold off; ylim([0 900000]); xlabel('Distance from Urban Area')   % Distance from Urban Area

end

%% Urbanization of soil/sand data
function [] = soil_sand( scene_data, scene_data_size )
    urban_075 = (scene_data(:,1) == 0.75 & scene_data(:,2) == 0.25);
    non_urban_075 = (scene_data(:,1) == 0.75 & scene_data(:,2) ~= 0.25);

    urban_075 = repmat(urban_075,1,scene_data_size(2));
    non_urban_075 = repmat(non_urban_075,1,scene_data_size(2));

    urban_075 = urban_075 .* scene_data;
    urban_075(all(urban_075==0,2),:) = [];
    non_urban_075 = non_urban_075 .* scene_data;
    non_urban_075(all(non_urban_075==0,2),:) = [];

    num_bins = 10;

    figure('name','Urbanization of soil/sand (0.75) data');
    subplot(1,7,1); histogram(scene_data(:,1),4); hold on; histogram(scene_data(:,2),4); hold off; xlabel('Land Useage')                                % Land Use
    legend('urban','non-urban','Location','northwest');
    subplot(1,7,2); histogram(urban_075(:,3),num_bins); hold on; histogram(non_urban_075(:,3),num_bins); hold off; xlabel('Slope')                      % Slope
    subplot(1,7,3); histogram(urban_075(:,4),num_bins); hold on; histogram(non_urban_075(:,4),num_bins); hold off; xlabel('Distance from I-25 & I-40')	% Distance from I25/I40
    subplot(1,7,4); histogram(urban_075(:,5),num_bins); hold on; histogram(non_urban_075(:,5),num_bins); hold off; xlabel('Distance from Rio Grande')	% Distance from Rio Grande
    subplot(1,7,5); histogram(urban_075(:,6),num_bins); hold on; histogram(non_urban_075(:,6),num_bins); hold off; xlabel('Distance from City Center')	% Distance from City Center
    subplot(1,7,6); histogram(urban_075(:,7),num_bins); hold on; histogram(non_urban_075(:,7),num_bins); hold off; xlabel('Distance from Roads')        % Distance from Roads
    subplot(1,7,7); histogram(urban_075(:,8),num_bins); hold on; histogram(non_urban_075(:,8),num_bins); hold off; ylim([0 900000]); xlabel('Distance from Urban Area')   % Distance from Roads
end
