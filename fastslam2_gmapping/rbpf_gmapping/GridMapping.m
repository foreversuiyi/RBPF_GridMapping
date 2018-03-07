function [P_new, L_new] = GridMapping(L_past, pose, scan, UsableArea, L0, gridsize, LRF) 
%% Calculation of the grid map
    % Find area of interest
    x = round(pose(1));
    y = round(pose(2));
    L_new = L_past;
    
    for a = (x - UsableArea) : (x + UsableArea)
         if a >= 1 && a <= gridsize(1)
            for b = (y - UsableArea) : (y + UsableArea)
                if b >= 1 && b <= gridsize(2)
                    L_new(gridsize(2)-b, a) = L_past(gridsize(2)-b, a) + ...
                        InverseLaserModel(a,b,pose,scan,L0,UsableArea,LRF) - L0;
                end
            end
         end
    end
    % Calculate normalized probablity form log-odds
    P_new = 1 - (1./(1 + exp(L_new)));
end
