%% CA Model
function [ sim_img ] = ca( pdev, scene, nhood_size )

    A = 0.2;                   % Adjusting factor 
    scene_size = size(scene);  % Determine the size of the scene
    
    % Perform the celluar automata (CA) 
    % pdev = prob_dev(scene_Ag, Abm, nhood_size, k);
    con = combined_constraint(scene);
    omega = nlfilter((scene == 0.25),[nhood_size nhood_size],@omega_fun);
    pt = A .* reshape(pdev,scene_size(1),scene_size(2)) .* con .* omega;
    x = rand(scene_size(1),scene_size(2));
    urb = (pt > x) * 0.25;
    non_urb = (pt < x) .* scene;
    sim_img = urb + non_urb;
end

function [ con ] = combined_constraint( scene )
    % Constraints such as rivers, steep slopes, and protected lands can be 
    % incorporated into the probability estimation; defined as "con"
    con = scene ~= 1 | scene ~= 0;
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

