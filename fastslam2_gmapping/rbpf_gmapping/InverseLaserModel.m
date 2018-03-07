function result = InverseLaserModel(ii, jj, pose, scan, L0, UsableArea, LRF)
%% Calculation of the inverse model returning the likelihood of occupancy
    alpha= 3;               % Thickness of obsatcles
    beta= 0.2;              % Width of sensor beam
    % Center of mass of cell
    xc = ii-0.5;
    yc = jj-0.5;
    % Calculate the Distance Between the Cell and Robot
    r = sqrt((xc-pose(1))^2 + (yc-pose(2))^2);  
    theta = atan2(yc-pose(2), xc-pose(1)) - pose(3);
    theta = NormalizeAngle(theta);
    %Find estimate of number of scan closer to theta
    max_angle = LRF.MaxAngle * pi /180;
    rng=floor((theta + max_angle)/(2*max_angle/(LRF.FOV/LRF.AngleResolution +1)));
    
    if rng >= size(scan,2)-1
        rng = size(scan,2)-1;
    end
    if rng < 2
        rng = 2;
    end
    %Search around the number of scan
    argmin = abs(theta - scan(1,1));
    kmin=1;
    for d = rng-1:rng+1
       a = abs(theta - scan(1,d));
       if a < argmin
           argmin = a;
           kmin = d;
       end
    end
    if (r > min(UsableArea, scan(2,kmin)+alpha/2)) || (abs(theta-scan(1,kmin)) > beta/2)
        result = L0;
    elseif (scan(2,kmin) < UsableArea) && (abs(r - scan(2,kmin)) < alpha/2)
        result = 1.5;
    elseif (r <= scan(2,kmin))
        result = -1.5;
    else
        error('Grid occupancy failed');
    end
end