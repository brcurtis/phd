clear variables;
clc;

% Default input path for the imagery
inputPath = '/Users/briancurtis/Documents/PhD/model/data/imagery/';

% Read imagery
imagery = ['1990133';'1995163';'2016157'];
img01_date = '1990133';
img02_date = '2016157';
filename_extent = '_classified_aoi.tif';
img01 = strcat(inputPath,img01_date,'/',img01_date,filename_extent);
img02 = strcat(inputPath,img02_date,'/',img02_date,filename_extent);
A = imread(img01,'TIF');
B = imread(img02,'TIF');

A_urban = A == 1;
B_urban = B == 1;

% Display images
figure
subplot(1,2,1); imshow(A_urban,[]); title(img01_date);
subplot(1,2,2); imshow(B_urban,[]); title(img02_date);

diff_urban = sum(B_urban(:)) / sum(A_urban(:));