function [ ] = gen_gis_plots( gis, gis_data, inputPath_gis )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    for i = 1:gis_size(1)
        figure;
        colormap('jet');
        imagesc(gis_data(:,:,i));
        %colorbar;
        set(gca,'FontSize',24,'xtick',[],'ytick',[]);
        title(gis(i));
    end
    

end

