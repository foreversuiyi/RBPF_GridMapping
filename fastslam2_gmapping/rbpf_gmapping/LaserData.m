function scan = LaserData(pose, map, LRF)
% Extract Sensor Parameters
max_distance        = LRF.MaxDistance;
angle_resolution    = LRF.AngleResolution;
fov                 = LRF.FOV;
sigma               = LRF.Sigma;

% Extract Pose
x = pose(1);
y = pose(2);
theta = pose(3);

% Number of Scanning Lines
num_scans = floor(fov/angle_resolution) + 1;
rad_resolution = angle_resolution * pi / 180.0;
num_walls = size(map, 2);

% Global System Coordinates
x_start = zeros(1, num_walls);
y_start = zeros(1, num_walls);
x_end = zeros(1, num_walls);
y_end = zeros(1, num_walls);

% Laser Scanner Local System Coordinates
trans_x_start = zeros(1, num_walls);
trans_y_start = zeros(1, num_walls);
trans_x_end = zeros(1, num_walls);
trans_y_end = zeros(1, num_walls);

a = zeros(1, num_walls);
b = zeros(1, num_walls);
c = zeros(1, num_walls);
scan = zeros(2, num_walls);

% Start and End Values for the Lines Ending Points
for i = 1:num_walls
    x_start(i) = map(1,i);
    y_start(i) = map(2,i);
    x_end(i)   = map(3,i);
    y_end(i)   = map(4,i);
% Transform the lines to the laser scanner's coordinate system.
%[trans_x]   [ cos sin  -cos*x-sin*y ] [orig_x]
%[trans_y] = [-sin cos  sin*x -cos*y ] [orig_y]
%[      1]   [   0   0             1 ] [     1] 
    trans_x_start(i)=(x_start(i)-x)*cos(theta) + (y_start(i)-y)*sin(theta);
    trans_y_start(i)=(y_start(i)-y)*cos(theta) - (x_start(i)-x)*sin(theta); 
    trans_x_end(i) = (x_end(i)-x)*cos(theta) + (y_end(i)-y)*sin(theta);
    trans_y_end(i) = (y_end(i)-y)*cos(theta) - (x_end(i)-x)*sin(theta);
    % Starting and ending points are swapped if the x value of the 
    % starting point is bigger than the x value of the ending point.
    if trans_x_start(i) > trans_x_end(i)
        trans_x_temp = trans_x_start(i);
        trans_x_start(i) = trans_x_end(i);
        trans_x_end(i) = trans_x_temp;
        trans_y_temp = trans_y_start(i);
        trans_y_start(i) = trans_y_end(i);
        trans_y_end(i) = trans_y_temp;
    end;
% The line equations a*x+ b*y=c
  a(i)=trans_y_start(i)-trans_y_end(i);
  b(i)=trans_x_end(i)-trans_x_start(i);
  l=sqrt(a(i)^2+b(i)^2);
  a(i)=a(i)/l;
  b(i)=b(i)/l;
  c(i)=a(i)*trans_x_start(i)+b(i)*trans_y_start(i);
end;
% Convert to polar coordinates
% For each laser scanner angle
for i = 1:num_scans
    % Laser scanner maximum measured distance maxDistance(in meters)
    % Closest distance from the laser scanner.
    min_dist = max_distance;
    % Current laser scanner angle 
    phi = - fov*pi/180.0/2.0 + (i-1)*rad_resolution;   
    % Find the minimum distance from the lines to laser scanner in the current angle
    cos_phi = cos(phi);
    sin_phi = sin(phi);
    for j = 1:num_walls 
     temp = a(j)*cos_phi + b(j)*sin_phi;
     if ( abs(temp) > 1e-6)
        t = c(j)/temp;  % Signed length from 0 to the intersection, x/cos_phi
        if (t > 0 && t < min_dist)
           if (abs(trans_x_start(j) - trans_x_end(j))>1e-6 )
               if(t*cos_phi < trans_x_end(j) && t*cos_phi > trans_x_start(j))
                 min_dist=t;
               end
           else
             if (trans_y_end(j) > trans_y_start(j) )
                if(t*sin_phi < trans_y_end(j) && t*sin_phi > trans_y_start(j))
                  min_dist=t; 
                end
             else
                if(t*sin_phi > trans_y_end(j) && t*sin_phi < trans_y_start(j))
                  min_dist=t; 
                end
             end
           end
        end
      end    
    end
    % Return the Simulated Laser Scan Data
    scan(1,i) = NormalizeAngle(phi + sigma(1)*randn);
    scan(2,i) = min_dist + sigma(2)*randn;
end