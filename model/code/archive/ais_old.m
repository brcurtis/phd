function [ Ab, Af_final ] = ais( scene_size, Ag_num_samples, Ab_init_size, Ngen, scenes, slope, distancefrom_I25_I40, distancefrom_rio_grande, distancefrom_city_center )
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

    %% Generate the inital antigen library
    % Create the random sample point matrix (Ag_samples) containing the 
    % x and y coordinates
    Ag_samples(:,1) = randi([1 scene_size(1)],round((scene_size(1) * scene_size(2) * Ag_num_samples)),1);
    Ag_samples(:,2) = randi([1 scene_size(2)],round((scene_size(1) * scene_size(2) * Ag_num_samples)),1);

    size_Ag = size(Ag_samples);

    for i=1:size_Ag(1)
        Ag(i,1) = scenes(Ag_samples(i,1),Ag_samples(i,2),1);
        Ag(i,2) = slope(Ag_samples(i,1),Ag_samples(i,2));
        Ag(i,3) = distancefrom_I25_I40(Ag_samples(i,1),Ag_samples(i,2));
        Ag(i,4) = distancefrom_rio_grande(Ag_samples(i,1),Ag_samples(i,2));
        Ag(i,5) = distancefrom_city_center(Ag_samples(i,1),Ag_samples(i,2));
        
        % Urbanization results from comparing the first scene to the
        % second scene (i.e. did the sample go from non-urban to urban?)
        if scenes(Ag_samples(i,1),Ag_samples(i,2),1) ~= 0.25 && scenes(Ag_samples(i,1),Ag_samples(i,2),2) == 0.25
            Ag(i,6) = 1;
        else
            Ag(i,6) = 0;
        end 
    end
    
    % Determine the library size for later use
    Ag_size = size(Ag);
    
    %% Randomly genterate the initial antibody library (Ab)
    %Ab = generate(round((scene_size(1) * scene_size(2) * Ab_init_size)));
    Ab = generate(Ab_init_size);
    Ab_size = size(Ab);
    
    %% AIS 
    % [--- insert text here ---]
    % Setup
    ttime = 0;
    
    % Inputs
    n = round(Ab_size(1) * 0.50);       % n number of antibodies to select and clone with highest affinities
    beta = 0.50;                        % multiplying factor for cloning
    b = 1.5;                              % constant for controlling antibody mutation probability
    d = 0.50;                           % constant for replacing the d(%) lowest antibodies in Ab with new randomly generated antibodies
    
    % Perform the AIS-learing by looping through a user-defined number of generations
    for t = 1:Ngen
        tic
        for i = 1:Ag_size(1)
            Af          = affinity(Ab, Ag(i,:));                % Step 1: Measure affinity (Af) of the antigen (Ag(i,:)) against each antibody (Ab)
            Abn         = select(Ab, Af, n);                    % Step 2: Select top n antibodies (Ab)
            C           = clone(Abn, beta, n);                  % Step 3: Clone the top n antibodies (Abn)
            C_mut       = mutation(C, Ag(i,:), b);              % Step 4: Mutate newly cloned anitbodies (C)
            Af_mut      = affinity(C_mut, Ag(i,:));             % Step 5: Mesaure affinity of the new anitbodies against the antigen (Ag(i,:)) 
            Abn         = select(C_mut, Af_mut, n);             % Step 6: Select the highest ranking antibodies
            Ab          = insert(Ab, Af, Abn);                  % Step 7: Insert result into antibody library
            temp(i,:)   = max(Ab(:,7));
            Abd         = generate(round(d*Ab_size(1)));        % Step 8: Generate d new anitbodies to replace the lowest affinity antibodies
            Ab          = replace(Ab, Abd);                     % Step 9: Replace the d lowest affinities in Ab with new antibodies (Abn)            
        end
        Af_final(:,t) = temp;
        fprintf('Completed Ngen %d in %f seconds\n',t,toc);
        ttime = ttime + toc;
    end
    fprintf('Total time: %f seconds\n',ttime);
end


%% Support functions

function [ Af ] = affinity( Ab, Ag );
    Ab_size = size(Ab);
    for i = 1:Ab_size(1)
        % Determine the fitness    
        Af(i) = 1 / (1 + (sqrt(sum((Ag(:,1:Ab_size(2)-1) - Ab(i,1:Ab_size(2)-1)).^2))));
        %Af(i) = 1 / (1 + (sqrt(sum((Ag - Ab(i,1:Ab_size(2))).^2))));
    end
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
    C_mut(:,1) = C(:,1);                        % Copy in the first attribute (land use); not mutated
    C_mut(:,C_size(2)-1) = C(:,C_size(2)-1);    % Copy in the last attribute (conversion results); not mutated
end

% Insert results into the antibodies library (Ab)
function [ Ab ] = insert( Ab, Af, Abn )
    Ab_org_size = size(Ab);
    Ab(:,Ab_org_size(2)+1) = Af;            % Append affinity onto Ab
    Ab = [Ab;Abn];                          % Concatenate Abn into Ab
    Ab_size = size(Ab);                     
    Ab = sortrows(Ab,-Ab_size(2));          % Order anitbodies (Ab) by their affinity (Af)
    Ab(Ab_org_size(1) + 1:end,:) = [];      % Truncate Ab to original size (Ab_org_size)
end

% Generate d antibodies for (1) initialization and (2) for new antibodies 
% to replace the lowest affinity antibodies in Ab
function [ Abd ] = generate( d )
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






%% Throw-away code
%             % Determine the level of concentration (Lcon)
%             threshold = 0.8; % threshold value for the affinity
%             temp = Af >= threshold;
%             Lcon = sum(temp,1);
%             %x = size(Lcon);
%             %Lcon = repmat(Lcon,1,x(1));
% 
%             % Determine the cloning rate (Pg)
%             alpha = 2.0;
%             Pg = (alpha .* Af) ./ Lcon;
%             Pg(isinf(Pg)) = 0;
% 
%             % Mutation
%             %for j = 1:Ab_lib_size(1)
%             %    Ab_lib(j,:) - ((1-exp(sqrt(sum((Ag_lib(1,:) - Ab_lib(1,:)).^2)))) * (Ab_lib(1,:) - Ag_lib(1,:)));



%     Ab(:,1) = randi([1 4],round((scene_size(1) * scene_size(2) * Ab_init_size)),1) * .25;
%     Ab(:,2) = rand(round((scene_size(1) * scene_size(2) * Ab_init_size)),1);
%     Ab(:,3) = rand(round((scene_size(1) * scene_size(2) * Ab_init_size)),1);
%     Ab(:,4) = rand(round((scene_size(1) * scene_size(2) * Ab_init_size)),1);
%     Ab(:,5) = rand(round((scene_size(1) * scene_size(2) * Ab_init_size)),1);
%     Ab(:,6) = randi([0 1],round((scene_size(1) * scene_size(2) * Ab_init_size)),1);