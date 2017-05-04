%% Setup
clear variables;
clc;
addpath('/Users/briancurtis/Documents/PhD/model/code');

%% Load data
% Locations and gis variable for Albuquerque, NM
inputPath_gis = '/Users/briancurtis/Documents/PhD/model/data/gis_data/alb/';
inputPath_img = '/Users/briancurtis/Documents/PhD/model/data/imagery/alb/';
images = [ '2001';'2006';'2011' ];
gis = { 'slope';'expressways';'major_water';'city_centers';'roads';'parks';'schools';'population';'railways'};

% Locations for Anderson, Pickens, and Greenville County, SC
% inputPath_gis = '/Users/briancurtis/Documents/PhD/model/data/gis_data/sc/';
% inputPath_img = '/Users/briancurtis/Documents/PhD/model/data/imagery/sc/';
% images = [ '2001';'2006';'2011' ];
% gis = { 'slope';'expressways';'major_water';'city_centers';'roads';'parks';'schools';'population';'railways' };

% Locations for Jackson County, OR
% inputPath_gis = '/Users/briancurtis/Documents/PhD/model/data/gis_data/oregon/';
% inputPath_img = '/Users/briancurtis/Documents/PhD/model/data/imagery/oregon/';
% images = [ '2001';'2006';'2011' ];
% gis = { 'slope';'expressways';'major_water';'city_centers';'roads';'parks';'schools';'population';'railways' };

image_combo = [1,2;1,3;2,3];

[ imagery, gis_data ] = load_data( images, inputPath_img, gis, inputPath_gis );

%% Setup variables and sample set generation
num_samples = 0.004;    % Percentage of image size to sample
[ all_data, sample_set ] = genSamples( imagery, image_combo, gis_data, num_samples );   % Generate the samples and prepare the data

%% Machine Learning (ANN) - Trainer
predictors = sample_set(:,[1,4:13])';
response = sample_set(:,[14])';
trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation

% Create a Pattern Recognition Network
hiddenLayerSize = 10;
net = patternnet(hiddenLayerSize, trainFcn);
net.trainParam.showWindow = false;              % hide the training window

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,predictors,response);


%% The Model
nhood_size = 3;               % Neighborhood size for CA (Omega)

% Identify urban pixels in each image (2001, 2006, and 2011)
img_2001 = imagery(:,:,1) >= 22 & imagery(:,:,1) <= 24;
img_2006 = imagery(:,:,2) >= 22 & imagery(:,:,2) <= 24;
img_2011 = imagery(:,:,3) >= 22 & imagery(:,:,3) <= 24;

% Setup the comparison image for FoM metric
comp_img = img_2001 ~=1 & img_2006 == 1;

% Starting with 2001 imagery
curr_img = imagery(:,:,1);
curr_img = curr_img >= 22 & curr_img <= 24;
size_curr_img = size(curr_img);

% Define number of elements to keep per loop
% num_elem = 100;
num_elem_perc = 0.005;
num_elem = ((sum(sum((img_2006 == 1))) - sum(sum((img_2001 == 1)))));
num_elem = round(num_elem * num_elem_perc);


% Define the number of loops based on the growth between 2006 and 2001
% num_loops = ((sum(sum((img_2006 == 1))) - sum(sum((img_2001 == 1)))) / num_elem);
% num_loops = round(num_loops + (num_loops * .20));

% Set up counter for the looping
FoM_end_counter = 0;
FoM_max = 0;
FoM_max_end = 10;
i = 1;

%for i = 1:num_loops
while FoM_end_counter ~= FoM_max_end
    tic
    % Set up the variables
    omega_urb = nlfilter(curr_img,[nhood_size nhood_size],@omega_fun);
    omega_urb = reshape(omega_urb,(size_curr_img(1)*size_curr_img(2)),1) .* 8;
    
    % Set up to run through machine learning algorithm
    curr_img = reshape(curr_img,(size_curr_img(1) * size_curr_img(2)),1);
    x = cat(2,curr_img,all_data(:,[4:12]),omega_urb)';
    
    % Run through machine learning algorithm and get results (probability)
    temp = net(x);
    
    % Sort and select the top 'n' non-urbanized cells from the network results
    temp = temp .* (curr_img ~= 1)';
    [temp_sort,temp_ind] = sort(temp,'descend');
    temp_ind = temp_ind(temp_ind ~= 0);
    temp_ind(num_elem+1:end) = [];
    curr_img(temp_ind) = 1;
    
    % Collect metrics
    %growth(i) = sum(curr_img) / sum(sum(imagery(:,:,1) >= 22 & imagery(:,:,1) <= 24));
    
    temp = [];
    temp = reshape(curr_img,size_curr_img(1),size_curr_img(2));
    temp = temp == 1 & img_2001 ~=1;
    % FoM(i) = figureofmerit( imagery(:,:,2), reshape(curr_img,size_curr_img(1),size_curr_img(2)));
    % FoM(i) = figureofmerit( comp_img, temp);
    FoM = figureofmerit( comp_img, temp);
    % c_matrix(:,:,i) = confusionmat(reshape(img_2006,(size_curr_img(1) * size_curr_img(2)),1),curr_img);
    c_matrix = confusionmat(reshape(img_2006,(size_curr_img(1) * size_curr_img(2)),1),curr_img);
    
    % Reshape curr_img for next pass through the loop's nlfilter (i.e.
    % needs to be 2-D matrix)
    curr_img = reshape(curr_img,size_curr_img(1),size_curr_img(2));
    % sim_img_store(:,:,i) = curr_img;
    
    if i ~= 1;
        if FoM < FoM_max
            FoM_end_counter = FoM_end_counter + 1;
            fprintf('** Completed %d in %f seconds with an FoM of %f **\n',i,toc,FoM);
        else
            % Store best simulated image
            sim_img_store = curr_img;
            
            % Store metrics
            FoM_max = FoM;
            FoM_end_counter = 0;
            c_matrix_max = c_matrix;
            
            %fprintf('Completed %d of %d in %f seconds with an FoM of %f\n',i,num_loops,toc,FoM(i));
            fprintf('Completed %d in %f seconds with an FoM of %f\n',i,toc,FoM);
        end
    end
    
    i = i + 1;
    
end

%accuracy = reshape((c_matrix(1,1,:) + c_matrix(2,2,i)) ./ (sum(sum(c_matrix(:,:,:)))),num_loops,1,1);

%% Code Wasteland
%images = [ '1990133';'1995163';'2016157' ];
%gis = { 'slope';'expressways';'major_water';'city_centers';'roads';'urban_area';'parks';'schools';'cultural_facilities';'population';'railways';'farmland'; 'industrial' };

% for i = 1:num_loops
%     tic
%     % Machine Learning - Predictor
%     sim_data = all_data;
%     sim_data(:,1) = reshape(curr_img,(size_curr_img(1) * size_curr_img(2)),1);
%     [ yfit, score ] = trainedClassifier.predictFcn(sim_data(:,[1,4:16]));
%     pdev(:,i) = score(:,2);
%     
%     % Cellular Automata
%     [ sim_img ] = ca( pdev(:,i), curr_img, nhood_size );
%     
%     % Collection metrics
%     growth(i) = sum(sum(sim_img == 0.25)) / sum(sum(imagery(:,:,1) == 0.25));
%     FoM(i) = figureofmerit( imagery(:,:,2), sim_img );
%     c_matrix(:,:,i) = confusionmat(reshape(imagery(:,:,2) == 0.25,(size_curr_img(1) * size_curr_img(2)),1),reshape(sim_img == 0.25,(size_curr_img(1) * size_curr_img(2)),1));
%     
%     % Set up for next pass through the loop
%     curr_img = sim_img;
%     sim_img_store(:,:,i) = sim_img;
%     
%     fprintf('Completed %d of %d in %f seconds with an FoM of %f\n',i,num_loops,toc,FoM(i));
% end