%% Clean
clear
close all
clc
addpath(genpath('rbpf_gmapping'))
rng(10)
%% Setup Parameters
% General Parameters
ts                  = 1;                    % Sampling Time             (unit: s)
num_cells           = 10;                   % Number of Cells           (unit: number per meter)
num_particles       = 3;                    % Number of Particles       (unit: number)
num_samples         = 10;                    % Number of Samples         (unit: number)
cov_odom            = [0.01 0; 0 0.001];    % Odometry Covarience       (unit: velocity)
max_speed           = 4;                    % Max_speed * num_cells     (unit: cells/timestep)
usable_area         = 4*num_cells;          % Radius Usable Data        (unit: cells)
% Controller Parameters
con_kp_v            = 0.4;                  % Velocity Control Command
con_kp_w            = 0.8;                  % Angular Velocity Control Command
% Motion Model Velocity Error Parameters  
a=[0.1 0.001 0.001 0.1 0.0001 0.0001]; 
% Laser Range Finder
LRF.MaxDistance     = 20*num_cells;         % Maximum Measurement Distance  (unit: cells)
LRF.AngleResolution = 0.5;                  % Measurement Angle Resolution  (unit: degree)
LRF.FOV             = 180;                  % Measurement Field of View     (unit: degree)
LRF.MaxAngle        = 90;                   % Maximum Measurement Angle     (unit: degree)
LRF.Sigma           = [0.01, 0.01];         % Range Sensor Deviation        (unit: degree; cells)

%% Generate Map
% Exterior Map
exterior_x = [150,225,225,250,250,150,150,25,25,100,100,150];
exterior_y = [25,25,75,75,200,200,275,275,75,75,25,25];
% Interior Map
interior_x_1 = [125,200,200,125,125];
interior_y_1 = [50,50,75,75,50];
interior_x_2 = [50,125,125,50,50];
interior_y_2 = [200,200,250,250,200];
interior_x_3 = [50,100,100,50,50];
interior_y_3 = [100,100,175,175,100];
interior_x_4 = [125,225,225,125,125];
interior_y_4 = [100,100,175,175,100];
map = [ exterior_x(1:end-1),interior_x_1(1:end-1),interior_x_2(1:end-1),interior_x_3(1:end-1),interior_x_4(1:end-1);...
        exterior_y(1:end-1),interior_y_1(1:end-1),interior_y_2(1:end-1),interior_y_3(1:end-1),interior_y_4(1:end-1);...
        exterior_x(2:end),interior_x_1(2:end),interior_x_2(2:end),interior_x_3(2:end),interior_x_4(2:end);...
        exterior_y(2:end),interior_y_1(2:end),interior_y_2(2:end),interior_y_3(2:end),interior_y_4(2:end)];

via_points = [110,40;200,40;210,50;220,90;240,100;240,180;220,190;150,190;140,200;...
    140,250;130,260;50,260;40,250;40,200;50,190;100,190;110,180;110,100;100,90;50,90;...
    40,100;40,180;50,190;220,190;240,180;240,100;230,90;120,90;110,80;110,40];
% Generate Path Using Points and Desired Velocity as Boundary Conditions
path = mstraj(via_points,[max_speed, max_speed],[],[via_points(1,1) via_points(1,2)],ts,0);

%% Grid Occupancy Initialization
num_steps   = size(path, 1);
bound_scale = 1.2;                                      % The Scale of Boundary to The Map
map_length  = max(max(map)) * bound_scale;              % Length of the Grid
map_width   = map_length;                               % Width of the Grid 
init_log    = zeros(map_width, map_length);             % Initialize the Log Probability of the Map
%% Data Initialization
% Velocity commands and robot states
con_cmd = zeros(2, num_steps);
x_states = zeros(3, num_steps);
x0 = [via_points(1,1), via_points(1,2), 0];
% Particle Filter Structure
PF.Iter                 = 1;                        % Iteration Number
PF.NumP                 = num_particles;            % Number of Particles
PF.Weight               = zeros(PF.NumP, num_steps);% Particle Weights
PF.Sigma                = LRF.Sigma;                % The noise of Sensor
PF.Particles            = zeros(3, PF.NumP, num_steps);             % Particles Position
PF.Probability          = zeros(map_width, map_length, PF.NumP);    % Probability of Grid Map
PF.Log                  = zeros(map_width, map_length, PF.NumP);  % Log Odds Probability
PF.Log0                 = 0;                        % L0 for Inverse Range Sensor Model
PF.GridSize             = [map_length, map_width,]; % Grid Size
PF.Xstate               = zeros(3, num_steps);      % Estimate Position
PF.CovOdom              = cov_odom;                 % Variance Matrix of Odometry
PF.UsableArea           = usable_area;              % Laser Scanner Usable Area
PF.Best                 = ones(1, num_steps);       % Best Particle Index with Highest Weight
%% Initialization
tic;
x_states(:, 1) = x0;                                    % Initialize Particle Position
PF.Xstate(:, 1) = x0;                                  % Initialize Estimate
laser_data = LaserData(PF.Xstate(:, 1), map, LRF);
% Initialize the Grid Map
[~, L_initial] = GridMapping(init_log, PF.Xstate(:, 1), laser_data, ...
    PF.UsableArea, PF.Log0, PF.GridSize, LRF);
[P_initial, L_initial] = GridMapping(L_initial, PF.Xstate(:, 1), laser_data, ...
    PF.UsableArea, PF.Log0, PF.GridSize, LRF);
% Intialize Variables
for i = 1:num_particles
    PF.Probability(:,:,i) = P_initial;
    PF.Log(:,:,i) = L_initial;
    PF.Weight(i, 1) = 1/num_particles;
    PF.Particles(:, i, 1) = x0;
end

%% Plot Map and path
fig_mapping = figure;
set(gcf,'units','normalized','outerposition',[0.1 0.2 0.8 0.6]);
subplot(1,3,1)
hold on
for d = 1:size(map,2)
    plot([map(1,d) map(3,d)],[map(2,d) map(4,d)],'k')
end
grid on
plot(path(:, 1), path(:, 2), 'r')
axis([0 300 0 300])
axis square

%% Start Simulation
for k = 2:num_steps
    % 1. Compute Error Signals
    error = sqrt((path(k,1) - x_states(1,k-1))^2 + (path(k,2) - x_states(2,k-1))^2);
    th = atan2(path(k,2)-x_states(2,k-1), path(k,1)-x_states(1,k-1));
    error_th = th - x_states(3, k-1);
    error_th = mod(error_th + pi, 2*pi) - pi;
    % 2. Control Signal
    con_cmd(1,k) = con_kp_v*error;
    con_cmd(2,k) = con_kp_w*error_th;
    % 3. True Motion Model
    x_states(:,k) = MotionCommandModel(x_states(:,k-1),con_cmd(:,k),ts);
    % 4. LRF Data
    laser_data = LaserData(x_states(:,k),map,LRF);
    % 5. The Rao-Blackwellized Particle Filter
    PF.Iter = k;
    PF = RaoBlackPF(x_states(:,k),laser_data,LRF,con_cmd(:,k),a,PF,num_samples,ts,fig_mapping);
%% Plot the Mapping Procedure
    figure(fig_mapping)
    subplot(1,3,2)
    caxis([0, 1]);
    colormap gray;
    colormap(flipud(colormap));
    hold on;
    imagesc(flipud(PF.Probability(:,:,PF.Best(k))));
    plot(PF.Xstate(1,k), PF.Xstate(2,k), 'dk')
    plot(x_states(1,k), x_states(2,k), 'or')
    axis([0 300 0 300])
    axis square
    hold off
    pause(0.0001)
end
rmpath(genpath('rbpf_gmapping'))
toc
