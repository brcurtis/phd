%% CA Model
function [ scene_sim ] = ca( scene, Abm, slope, distancefrom_I25_I40, distancefrom_rio_grande, distancefrom_city_center )

    A = 1;                          % Adjusting factor
    nhood_size = 3;                 % Neighborhood size for CA (Omega)
    k = 20;                         % number of anitbodies to select for pdev
    
    scene_size = size(scene);
    
    % create a temp scene with only highlights urban (0), water (127) and
    % non-urban (255)
    x = 255 * (scene ~= 0.25 & scene ~= 1);
    y = 127 * (scene == 1);
    scene_temp = x + y;
    
    % Build the scene antigenlibrary as a 3-D matrix
    scene_Ag(:,:,1) = scene;
    scene_Ag(:,:,2) = slope;
    scene_Ag(:,:,3) = distancefrom_I25_I40;
    scene_Ag(:,:,4) = distancefrom_rio_grande;
    scene_Ag(:,:,5) = distancefrom_city_center;
    
    % Perform the celluar automata (CA) 
    pdev = prob_dev(scene_Ag, Abm, nhood_size, k);
    con = combined_constraint(scene);
    omega = nlfilter((scene_Ag(:,:,1) == 0.25),[nhood_size nhood_size],@omega_fun);
    pt = A .* pdev .* con .* omega;
    x = rand(scene_size(1),scene_size(2));
    urb = (pt > x) * 0.25;
    non_urb = (pt < x) .* scene_Ag(:,:,1);
    scene_sim = urb + non_urb;
end

function [ pdev ] = prob_dev( scene_Ag, Abm, nhood_size, k )   
    scene_Ag_size = size(scene_Ag);
    temp = reshape(scene_Ag,(scene_Ag_size(1) * scene_Ag_size(2)),scene_Ag_size(3));
    temp_size = size(temp);
    final_time = 0;
    tic
    for i=1:temp_size(1)
        if temp(i,1) ~= 0.25 && temp(i,1) ~= 0
            %Abm_size = size(Abm);
            Af = affinity(Abm, temp(i,:));                          % Determine the fitness (i.e. Affinity)
            Abn = select( Abm, Af, k );                             % Select the top K antibodies
            urb = sum((Abn(:,6) == 1) .* Abn(:,7));                 % Determine the sum of the affinities for urbanize anitbody cells
            non_urb = sum((Abn(:,6) == 0) .* Abn(:,7));             % Determine the sum of the affinities for non-urbanize anitbody cells
            pdev(i) = urb / (urb + non_urb);
        else
            pdev(i) = 0;
        end
        % Timing for debugging purposes
        if mod((i/scene_Ag_size(2)),1) == 0
            final_time = final_time + toc;
            fprintf('Row %d complete in %f seconds - %f%% complete \n',(i/scene_Ag_size(2)),toc,((i/temp_size(1))*100))
            tic
        end
    end
    pdev = reshape(pdev,scene_Ag_size(1),scene_Ag_size(2));
    fprintf('pdev completed in %f seconds \n',final_time)
end

function [ con ] = combined_constraint( scene )
    % Constraints such as rivers, steep slopes, and protected lands can be 
    % incorporated into the probability estimation; defined as "con"
    con = scene ~= 1;
end

function [ Af ] = affinity( Abm, scene_Ag )
    Abm_size = size(Abm);
    % Determine the fitness (i.e. Affinity)    
    x = repmat(scene_Ag,Abm_size(1),1);
    Af = 1 ./ (1 + sqrt(sum((x(:,1:Abm_size(2)-1) - Abm(:,1:Abm_size(2)-1)).^2,2)));
end

function [ Abn ] = select( Abm, Af, k)
    Abm_size = size(Abm);
    Abm(:,Abm_size(2)+1) = Af;          % Append affinity onto Ab
    Abm_size = size(Abm);
    Abm = sortrows(Abm,-Abm_size(2));	% Order anitbodies (Ab) by their affinity (Af)
    Abn = Abm(1:k,:);                   % Select the top n antibodies
end

function [ omega ] = omega_fun(scene)
    if sum(sum(scene)) == 0
        omega = 0;
    elseif scene(2,2) == 1
        omega = (sum(sum(scene)) - 1) / 8;
    else
        omega = (sum(sum(scene))) / 8;
    end
end

%% Staging code for departure


% function [ pdev ] = prob_dev( scene_Ag, Abm, nhood_size, k )   
%     final_time = 0;
%     scene_Ag_size = size(scene_Ag);
%     for i=((nhood_size+1)/2):scene_Ag_size(1)
%         tic
%         for j=((nhood_size+1)/2):scene_Ag_size(2)
%             if scene_Ag(i,j,1) ~= 0.25 && scene_Ag(i,j,1) ~= 0
%                 Abm_size = size(Abm);
%                 temp = reshape(scene_Ag(i,j,:),[1,scene_Ag_size(3)]);   % Reshape the scene_Ag matrix for use in Af function
%                 Af = affinity(Abm, temp);                               % Determine the fitness (i.e. Affinity)
%                 Abn = select( Abm, Af, k );                             % Select the top K antibodies
%                 urb = sum((Abn(:,6) == 1) .* Abn(:,7));                 % Determine the sum of the affinities for urbanize anitbody cells
%                 non_urb = sum((Abn(:,6) == 0) .* Abn(:,7));             % Determine the sum of the affinities for non-urbanize anitbody cells
%                 pdev(i,j) = urb / (urb + non_urb);
%             else
%                 pdev(i,j) = 0;
%             end
%         end
%         final_time = final_time + toc;
%         fprintf('Row %d complete in %f seconds - %f%% complete \n',i,toc,((i/scene_Ag_size(1))*100))
%     end
%     fprintf('pdev completed in %f seconds \n',final_time)  
% end