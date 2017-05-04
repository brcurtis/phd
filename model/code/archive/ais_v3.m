function [ Ab, Abm, Afm, eval ] = ais( scene_size, Ag_num_samples, Ab_init_size, Ngen, scenes, slope, distancefrom_I25_I40, distancefrom_rio_grande, distancefrom_city_center )
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
    Ag  = genAg( scene_size, Ag_num_samples, scenes, slope, distancefrom_I25_I40, distancefrom_rio_grande, distancefrom_city_center );      
    Ag_size = size(Ag); % Determine the library size for later use
    
    Ag_eval = genAg( scene_size, Ag_num_samples, scenes, slope, distancefrom_I25_I40, distancefrom_rio_grande, distancefrom_city_center );
    Ag_eval_size = size(Ag_eval); % Determine the library size for later use
    
    %% Randomly genterate the initial antibody library (Ab)
    %Ab = generate(round((scene_size(1) * scene_size(2) * Ab_init_size)));
    Ab = genAb(Ab_init_size);
    Ab_size = size(Ab);
    
    %% AIS 
    % [--- insert text here ---]
    % Setup
    ttime = 0;
    Abm = zeros(Ag_size(1),Ab_size(2)+1);
    
    % Inputs
    n = round(Ab_size(1) * 0.25);      % n number of antibodies to select and clone with highest affinities
    beta = 0.50;                       % multiplying factor for cloning
    b = 1.50;                          % constant for controlling antibody mutation probability
    d = 0.25;                          % constant for replacing the d(%) lowest antibodies in Ab with new randomly generated antibodies
    
    % Perform the AIS-learing by looping through a user-defined number of generations
    for t = 1:Ngen
        tic
        for i = 1:Ag_size(1)
            Af                  = affinity(Ab, Ag(i,:));                % Step 1: Measure affinity (Af) of the antigen (Ag(i,:)) against each antibody (Ab)
            Abn                 = select(Ab, Af, n);                    % Step 2: Select top n antibodies (Ab)
            C                   = clone(Abn, beta, n);                  % Step 3: Clone the top n antibodies (Abn)
            C_mut               = mutation(C, Ag(i,:), b);              % Step 4: Mutate newly cloned anitbodies (C)
            Af_mut              = affinity(C_mut, Ag(i,:));             % Step 5: Mesaure affinity of the new anitbodies against the antigen (Ag(i,:)) 
            Abn                 = select(C_mut, Af_mut, n);             % Step 6: Select the highest ranking antibodies
          [ Abm(i,:), Ab ]      = insert(Abm(i,:), Ab, Af, Abn);        % Step 7: Insert result into antibody library (Ab) and memory (Abm)
            Abd                 = genAb(round(d*Ab_size(1)));           % Step 8: Generate d new anitbodies to replace the lowest affinity antibodies
            Ab                  = replace(Ab, Abd);                     % Step 9: Replace the d lowest affinities in Ab with new antibodies (Abd)            
        end
        Afm(:,t) = Abm(:,end);                                          % Store the Ngen affinities for each antibody
        eval(t,:) = ais_eval( Ag_eval, Ag_eval_size, Abm(:,1:6) );      % Perform the validation of the antibody library
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
                C_mut(i,j) = (C(i,j) - ((1 - exp(-d)) * (C(i,j) - Ag(:,j)))); % * (1 - r);
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
function [ Abm, Ab ] = insert( Abm, Ab, Af, Abn )
    Ab_org_size = size(Ab);
    Ab(:,Ab_org_size(2)+1) = Af;            % Append affinity onto Ab
    Ab = [Ab;unique(Abn,'rows')];           % Concatenate unique rows from Abn into Ab
    Ab_size = size(Ab);                     
    Ab = sortrows(Ab,-Ab_size(2));          % Order anitbodies (Ab) by their affinity (Af)
    Ab(Ab_org_size(1) + 1:end,:) = [];      % Truncate Ab to original size (Ab_org_size)
    
    if Abm(7) < Ab(1,end)
        Abm = Ab(1,:);
    end
end

% Generate d antibodies for (1) initialization and (2) for new antibodies 
% to replace the lowest affinity antibodies in Ab
function [ Abd ] = genAb( d )
    Abd(:,1) = randi([1 4],d,1) * .25;
    Abd(:,2) = rand(d,1);
    Abd(:,3) = rand(d,1);
    Abd(:,4) = rand(d,1);
    Abd(:,5) = rand(d,1);
    Abd(:,6) = randi([0 1],d,1);
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

function [ Ag ] = genAg( scene_size, Ag_num_samples, scenes, slope, distancefrom_I25_I40, distancefrom_rio_grande, distancefrom_city_center )    
    sample_size = round((scene_size(1) * scene_size(2) * Ag_num_samples));
    urb = round(sample_size / 2);
    non_urb = sample_size - urb;
    
    i = 1;
    while urb > 0
        Ag_sample(1) = randi([1 scene_size(1)]);
        Ag_sample(2) = randi([1 scene_size(2)]);
        if scenes(Ag_sample(1),Ag_sample(2),1) ~= 0.25 && scenes(Ag_sample(1),Ag_sample(2),2) == 0.25
            Ag(i,1) = scenes(Ag_sample(1),Ag_sample(2),1);
            Ag(i,2) = slope(Ag_sample(1),Ag_sample(2));
            Ag(i,3) = distancefrom_I25_I40(Ag_sample(1),Ag_sample(2));
            Ag(i,4) = distancefrom_rio_grande(Ag_sample(1),Ag_sample(2));
            Ag(i,5) = distancefrom_city_center(Ag_sample(1),Ag_sample(2));
            Ag(i,6) = 1;
            urb = urb - 1;
            i = i + 1;
        end
    end
    
    i = size(Ag);
    i = i(1) + 1;
    
    while non_urb > 0
        Ag_sample(1) = randi([1 scene_size(1)]);
        Ag_sample(2) = randi([1 scene_size(2)]);
        if scenes(Ag_sample(1),Ag_sample(2),1) ~= 0.25 && scenes(Ag_sample(1),Ag_sample(2),2) ~= 0.25
            Ag(i,1) = scenes(Ag_sample(1),Ag_sample(2),1);
            Ag(i,2) = slope(Ag_sample(1),Ag_sample(2));
            Ag(i,3) = distancefrom_I25_I40(Ag_sample(1),Ag_sample(2));
            Ag(i,4) = distancefrom_rio_grande(Ag_sample(1),Ag_sample(2));
            Ag(i,5) = distancefrom_city_center(Ag_sample(1),Ag_sample(2));
            Ag(i,6) = 0;
            non_urb = non_urb - 1;
            i = i + 1;
        end
    end
    Ag = sortrows(Ag,6);
end

%% Performance validation function

function [ eval ] = ais_eval( Ag_eval, Ag_eval_size, Abm )
    k = 20;
    for i = 1:Ag_eval_size(1)
        Af                  = affinity(Abm, Ag_eval(i,:));      % Step 1: Measure affinity (Af) of the antigen (Ag_val(i,:)) against each antibody (Abm)
        Abn                 = select(Abm, Af, k);               % Step 2: Select top k antibody (Abm)
        if Ag_eval(i,6) == 1
            eval(i) = sum(Abn(:,6)) / k;                        % percent with correct urbanization as a '1'
        else
            eval(i) = (k - sum(Abn(:,6))) / k;                  % percent with correct urbanization as a '0'
        end
    end
end
