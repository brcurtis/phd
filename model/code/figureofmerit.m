function [ FoM ] = figureofmerit( real_img, sim_img )
% For cell-level aspect, the indicator of ?figure of merit? (FoM) 
% has been commonly used by many researchers. The indicator of ?FoM? is 
% actually a ratio, where the numerator is the number of instances that 
% are observed developed and correctly simulated as developed, while the 
% denominator is the total num- ber of instances excluding persistently 
% non-changed instances.

% A represents the error due to observed developed and simulated as persistence
% A = sum(sum((real_img >= 22 & real_img <= 24) & sim_img ~= 1));
A = sum(sum(real_img == 1 & sim_img ~= 1));

% B is the agreement due to observed developed and simulated as developed
% B = sum(sum((real_img >= 22 & real_img <= 24) & sim_img == 1));
B = sum(sum(real_img == 1 & sim_img == 1));

% C is defined as the error due to observed developed and simulated as incorrect
% gaining category.  As the CA models only simulate the change of states 
% from non-urban to urban, the value of C should be equal to 0
C = 0;

% D is the error due to observed persistence and simulated as developed
% D = sum(sum((real_img <= 22 & real_img >= 24) & sim_img == 1));
D = sum(sum(real_img ~= 1 & sim_img == 1));

% Calculate the overal Figure of Merit (FoM)
FoM = (B / (A+B+C+D));

end

