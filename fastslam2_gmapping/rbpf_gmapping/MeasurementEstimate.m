function scan_estimate = MeasurementEstimate(pose,probability_map,resolution,FOV, max_range, UsableArea)
%% Extracts laser scan data from an occupancy grid map
    probability_map=flipud(probability_map);          % Map must be flipped because of the coordinates in which it is created
    object_limit = 0.817;               % Probability that says that an object is in that cell very certainly
	
	% Extract Pose
	x=pose(1);
	y=pose(2);
    theta = NormalizeAngle(pose(3));
    
	% Extract Pose in grid coordinates
	x_grid = round(pose(1));
	y_grid = round(pose(2));
	grid_sizex=size(probability_map,1);
	grid_sizey=size(probability_map,2);

	% Number of scanning lines (deduced from scan width and resolution)
	num_of_scanning_lines=floor(FOV/resolution)+1;
	rad_resolution=resolution*pi/180.0;

    % Initialize scan estimate to maximum scan range
	scan_estimate=...
        [linspace(-rad_resolution*num_of_scanning_lines/2, ...
        rad_resolution*num_of_scanning_lines/2,...
        num_of_scanning_lines); ...
        max_range*ones(1, num_of_scanning_lines)]; % Make sure this matches with other scans

    for a = (x_grid - UsableArea):(x_grid + UsableArea)
        if a >= 1 && a <= grid_sizex
            for b = (y_grid - UsableArea):(y_grid + UsableArea)
                if b >= 1 && b <= grid_sizey
                    if probability_map(b,a) >= object_limit
                        %Center of mass of cell
                        xc = a-0.5;
                        yc = b-0.5;
                        theta_cell = atan2(yc-y,xc-x);
                        theta_scan = theta_cell - theta;
                        theta_scan = NormalizeAngle(theta_scan);

                        % Find estimate of number of scan closer to theta
                        rng=floor((theta_scan+rad_resolution*num_of_scanning_lines/2)/rad_resolution);

                        % Calculate distance
                        r=sqrt((xc-x_grid)^2+(yc-y_grid)^2);

                        % Calculates how many scans will go through a cell depending on how far it is
                        theta_rng=floor((atan(0.5/r)/rad_resolution)); 
                        
                        % Update that scan angle if its new distance is found to be samaller
                        for d=rng-theta_rng:rng+theta_rng
                            if d <= 0 || d >= size(scan_estimate,2)
                                %skip
                            elseif r < scan_estimate(2,d)
                                scan_estimate(2,d)=r;
                            end
                        end		
                    end					
                end
            end
        end
    end      
end