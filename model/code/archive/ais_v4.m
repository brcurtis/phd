function [ Abm, Afm_eval, eval ] = ais( Ag_num_samples, Ab_init_size, Ngen, imagery, image_combo, gis_data )
    %% AIS model setup

    % Build the Antigen library - expressed as a series of vectors (Ag(i)) 
    % consisting of two parts: attributes and state conversion result. 
    % Attributes includes variable like land use, proximity variables, slope, 
    % etc.  Below is the construct of this Antigen library vector
    % Ag(i) = [ land use type; slope; distancefrom_I25_I40; distancefrom_rio_grande; distancefrom_city_center; conversion_result ]
    % Note, the "conversion_result" is the boolean value which results from 
    % the analysis between the previous image (e.g. 1990) and the next image 
    % (e.g. 1995).  If land was converted, from non-urban to urban, it is 
    % receives a '1', else '0'.

    %% Generate the inital antigen library (Ag) and the evaluation library (Ag_eval) 
    Ag  = genAg( Ag_num_samples, imagery, image_combo, gis_data );      
    Ag_size = size(Ag); % Determine the library size for later use
    
    Ag_eval = genAg( scene_size, Ag_num_samples, scenes, slope, distancefrom_expressways, distancefrom_major_water, distancefrom_city_centers, distancefrom_roads );
    Ag_eval_size = size(Ag_eval); % Determine the library size for later use
    
    %% Randomly genterate the initial antibody library (Ab)
    %Ab = generate(round((scene_size(1) * scene_size(2) * Ab_init_size)));
    Ab = genAb(Ab_init_size);
    Ab_size = size(Ab);
    
    %% AIS 
    % [--- insert text here ---]
    % Setup
    ttime = 0;
    Abm_max_size = 400;
    Abm = genAb(50);
    Abm = zeros(1,Ab_size(2));
    Abm(:,end+1) = 0;
    
    % Inputs
    n = round(Ab_size(1) * 0.25);       % n number of antibodies to select and clone with highest affinities
    beta = 0.50;                        % multiplying factor for cloning
    b = 1.50;                           % constant for controlling antibody mutation probability
    d = 0.5;                            % constant for replacing the d(%) lowest antibodies in Ab with new randomly generated antibodies
    
    % Perform the AIS-learing by looping through a user-defined number of generations
    for t = 1:Ngen
        tic
        for i = 1:Ag_size(1)
            Ab = vertcat(Ab,Abm(:,1:7));
            Af                  = affinity(Ab, Ag(i,:));                            % Step 1: Measure affinity (Af) of the antigen (Ag(i,:)) against each antibody (Ab)
            Abn                 = select(Ab, Af, n);                                % Step 2: Select top n antibodies (Ab)
            C                   = clone(Abn, beta, n);                              % Step 3: Clone the top n antibodies (Abn)
            C_mut               = mutation(C, Ag(i,:), b);                          % Step 4: Mutate newly cloned anitbodies (C)
            Af_mut              = affinity(C_mut, Ag(i,:));                         % Step 5: Mesaure affinity of the new anitbodies against the antigen (Ag(i,:)) 
            Abn                 = select(C_mut, Af_mut, n);                         % Step 6: Select the highest ranking antibodies
          [ Abm, Ab ]           = insert(Abm, Abm_max_size, Ab, Af, Abn, Ag(i,:));  % Step 7: Insert result into antibody library (Ab) and memory (Abm)
            Abd                 = genAb(round(d*Ab_size(1)));                       % Step 8: Generate d new anitbodies to replace the lowest affinity antibodies
            Ab                  = replace(Ab, Abd);                                 % Step 9: Replace the d lowest affinities in Ab with new antibodies (Abd)            
        end
        Afm_eval(:,t) = Abm(:,end);                         % Store the Ngen affinities for each antibody
        eval(:,t) = ais_eval( Ag, Abm(:,1:7) );        % Perform the validation of the antibody library
        fprintf('Completed Ngen %d in %f seconds\n',t,toc);
        ttime = ttime + toc;
    end
    fprintf('Total time: %f seconds\n',ttime);
end


%% Support functions

function [ Af ] = affinity( Ab, Ag )
    Ab_size = size(Ab);
    % Determine the fitness (i.e. Affinity)    
    x = repmat(Ag,Ab_size(1),1);
    Af = 1 ./ (1 + sqrt(sum((x(:,1:Ab_size(2)-1) - Ab(:,1:Ab_size(2)-1)).^2,2)));      % without binary flag for urbanization
    %Af = 1 ./ (1 + sqrt(sum((x(:,1:Ab_size(2)) - Ab(:,1:Ab_size(2))).^2,2)));          % with binary flag for urbanization
end

% Select: Selects the top n antibodies and returns the actual antibody in order of 
% their affinity (Af)
function [ Abn ] = select( Ab, Af, n)
    Ab_size = size(Ab);
    Ab(:,Ab_size(2)+1) = Af;            % Append affinity onto Ab
    Ab_size = size(Ab);
    Ab = sortrows(Ab,-Ab_size(2));      % Order anitbodies (Ab) by their affinity (Af)
    Abn = Ab(1:n,:);                    % Select the top n antibodies
end

% Clone: 
function [ C ] = clone( Abn, beta, n )
    x = 0;
    for i = 1:n
        y = round((beta * n)/i);
        C(x+1:x+y,:) = repmat(Abn(i,:),[y,1]);
        x = x + y;
    end
end

% Mutation:
function [ C_mut ] = mutation( C, Ag, b )
    C_size = size(C);                           % Determine the size of C for looping purposes
    for i = 1:C_size(1)
        for j = 2:C_size(2)-2                   % ignore first attribute (land use) and last two (conversion results & affinity)
            r = rand();                         % Generate a random number
            if (C(i,end) / b) < r
                d = (1 / C(i,end)) - 1;
                C_mut(i,j) = (C(i,j) - ((1 - exp(-d)) * (C(i,j) - Ag(:,j)))) * (1 - r);
            else
                C_mut(i,j) = C(i,j);
            end
        end
    end
    %C_mut(:,1) = C(:,1);                        % Copy in the first attribute (land use); not mutated
    C_mut(:,1) = randi(4,C_size(1),1) * 0.25;
    C_mut(:,C_size(2)-1) = C(:,C_size(2)-1);    % Copy in the last attribute (conversion results); not mutated
    
    % Attempt to transform from nested loops to matrix math
    %r_matrix = rand(C_size(1),C_size(2)-3);
    %test = C(:,end) / b;
    %test = repmat(test,1,C_size(2)-3);
    %test = test > r_matrix;
    %d = (1 ./ C(:,end)) - 1;
    %d = repmat(d,1,C_size(2)-3);
    %C_temp = C(:,2:5);
    %Ag_temp = repmat(Ag,C_size(1),1);
    %Ag_temp(:,[1 6]) = [];
    %temp = C_temp - ((1 - exp(-d)) .* (C_temp - Ag_temp)) .* (1 - r_matrix);
    
end

% Insert results into the antibodies library (Ab) and memory (Abm)
function [ Abm, Ab ] = insert( Abm, Abm_max_size, Ab, Af, Abn, Ag )
    Ab_org_size = size(Ab);
    Abm_size = size(Abm);
    Ab_org_size(1) = Ab_org_size(1) - Abm_size(1);
    
    Ab(:,Ab_org_size(2)+1) = Af;            % Append affinity onto Ab
    %Abm(:,Ab_org_size(2)+1) = Afm;          % Append affinity onto Abm
    Ab = [Ab;unique(Abn,'rows')];           % Concatenate unique rows from Abn into Ab
    Ab_size = size(Ab);                     
    Ab = sortrows(Ab,-Ab_size(2));          % Order anitbodies (Ab) by their affinity (Af)
    Ab(Ab_org_size(1) + 1:end,:) = [];      % Truncate Ab to original size (Ab_org_size)
    
    % find the 
    temp = Ab(:,1) == Ag(:,1) & Ab(:,end-1) == Ag(:,end);
    temp = repmat(temp,1,Ab_size(2));
    temp = temp .* Ab;
    [max_af,pos] = max(temp);
    Abm(end+1,:) = max_af;
    
    if Abm_size(1) > Abm_max_size
        Abm = sortrows(Abm,-Abm_size(2));
        Abm(Abm_max_size + 1:end,:) = [];
    end
    
%     [Afm_max,pos] = max(Afm(:,end));       % Order anitbodies (Abm) by their affinity (Af)
%     if Afm(pos) < Ab(1,end) && Afm(pos) > Abm(pos,end)
%         Abm(pos,:) = Ab(1,:);
%     end
end

% Generate d antibodies for (1) initialization and (2) for new antibodies 
% to replace the lowest affinity antibodies in Ab
function [ Abd ] = genAb( d )
    Abd(:,1) = randi([1 4],d,1) * .25;
    Abd(:,2) = rand(d,1);
    Abd(:,3) = rand(d,1);
    Abd(:,4) = rand(d,1);
    Abd(:,5) = rand(d,1);
    Abd(:,6) = rand(d,1);
    Abd(:,7) = randi([0 1],d,1);
end

% Replace the d (%) lowest affinity antibodies in Ab with newly generated
% antibodies (Abn) from the gnerate function
function [ Ab ] = replace( Ab, Abd )
    
    % Order anitbodies (Ab) by their affinity (Af)
    Ab_size = size(Ab);
    Ab = sortrows(Ab,-Ab_size(2));
    
    % Replace the d (%) lowest with new (Abd)
    Abd_size = size(Abd);
    Ab_size = size(Ab);
    Ab(:,Ab_size(2)) = [];                                      % remove affinity measures
    Ab( [(Ab_size(1) - Abd_size(1) + 1) : end], : ) = Abd;
end

function [ Ag ] = genAg( Ag_num_samples, imagery, image_combo, gis_data )    
    
    %urb = round(sample_size / 2);
    %non_urb = sample_size - urb;
    
    imagery_size = size(imagery);
    gis_data_size = size(gis_data);
    sample_size = round((imagery_size(1) * imagery_size(2) * Ag_num_samples));
    
    for i = 1:imagery_size(3)
        temp(:,i) = reshape(imagery(:,:,i),(imagery_size(1) * imagery_size(2)),1);
        temp(:,i) = temp(:,i) == 0.25;
    end
    
    x = size(image_combo);
    for i = 1:x(1)
       imagery_delta(:,i) =  temp(:,image_combo(i,1)) ~= 1 & temp(:,image_combo(i,2)) == 1;
    end
    
    for i = 1:gis_data_size(3)
        gis_sample(:,i) = reshape(gis_data(:,:,i),(imagery_size(1) * imagery_size(2)),1);
    end
    
%     scene_data(:,1) = reshape(scenes(:,:,1),(scene_size(1) * scene_size(2)),1);
%     scene_data(:,2) = reshape(scenes(:,:,2),(scene_size(1) * scene_size(2)),1);
%     scene_data(:,3) = reshape(slope,(scene_size(1) * scene_size(2)),1);
%     scene_data(:,4) = reshape(distancefrom_I25_I40,(scene_size(1) * scene_size(2)),1);
%     scene_data(:,5) = reshape(distancefrom_rio_grande,(scene_size(1) * scene_size(2)),1);
%     scene_data(:,6) = reshape(distancefrom_city_center,(scene_size(1) * scene_size(2)),1);
%     scene_data(:,7) = reshape(distancefrom_roads,(scene_size(1) * scene_size(2)),1);
%     scene_data(:,8) = scene_data(:,1) ~= 0.25 & scene_data(:,2) == 0.25;
%     
%     scene_data_size = size(scene_data);
    
    land_use = unique(imagery(:,:));
    land_use_size = size(land_use);
    i = 1;
    
    for urb = 0:1
        for x = 3:land_use_size(1)-1
            temp = imagery(:,1) == land_use(x) & imagery_delta(:,1) == urb;
            temp = repmat(temp,1,scene_data_size(2)) .* scene_data;
            temp(all(temp==0,2),:)=[];
            temp_size = size(temp);
            rand_sample = randi([1 temp_size(1)],round(sample_size/4),1);
            temp = temp(rand_sample,:);
            if exist('Ag') == 1
                Ag = vertcat(Ag,temp);
            else
                Ag = temp;
            end
        end
    end
    
    Ag = Ag(:,[1 3:end]);   % remove the scene data from column 2
end

%% Performance validation function

function [ eval ] = ais_eval( Ag_eval, Abm )
    k = 20;
    Ag_eval_size = size(Ag_eval);
    for i = 1:Ag_eval_size(1)
        Af                  = affinity(Abm, Ag_eval(i,:));      % Step 1: Measure affinity (Af) of the antigen (Ag_eval(i,:)) against each antibody (Abm)
        Abn                 = select(Abm, Af, k);               % Step 2: Select top k antibody (Abm)
        if Ag_eval(i,end) == 1
            eval(i) = sum(Abn(:,end-1)) / k;                        % percent with correct urbanization as a '1'
        else
            eval(i) = (k - sum(Abn(:,end-1))) / k;                  % percent with correct urbanization as a '0'
        end
    end
end
