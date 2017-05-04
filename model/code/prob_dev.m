function [ pdev ] = prob_dev( scene_Ag, Abm, nhood_size, k )                                           % Determine scene size for looping
    if scene_Ag(i,j,1) ~= 0.25 && scene_Ag(i,j,1) ~= 0
        Abm_size = size(Abm);
        temp = reshape(scene_Ag(i,j,:),[1,scene_size(3)]);      % Reshape the scene_Ag matrix for use in Af function
        Af = affinity(Abm, temp);                               % Determine the fitness (i.e. Affinity)
        Abn = select( Abm, Af, k );                             % Select the top K antibodies
        urb = sum((Abn(:,6) == 1) .* Abn(:,7));                 % Determine the sum of the affinities for urbanize anitbody cells
        non_urb = sum((Abn(:,6) == 0) .* Abn(:,7));             % Determine the sum of the affinities for non-urbanize anitbody cells
        pdev(i,j) = urb / (urb + non_urb);
    else
        pdev = 0;
    end
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