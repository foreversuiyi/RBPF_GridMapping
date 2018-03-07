function PF = RaoBlackPF(x_true, laser_data, LRF, con_cmd, a, PF, num_samples, ts, fig_handle)
%% Rao-BlackWellization Particle Filter
    % x_true: True position of robot
    % laser_data: Laser scan data
    % LRF: Laser scan inforamtion
    % con_cmd: Movement control
    % a: Odometry parameters found via Motion_Model_Velocity_Test
    % PF: PF class with all particle filter information
    % num_samples: Number of samples for the Particle Filter
    % ts: Time step
%% Variable Initialization
    iter = PF.Iter;
    num_particles = PF.NumP;
    weight_past = PF.Weight(:,iter-1);
    x_past= PF.Particles(:,:,iter-1);
    % Initialize new variables
    x_initial_estimate=zeros(3, num_particles);           % Initial estimate from previous Xt-1 and Ut
    laser_estimate=zeros(2, 361, num_particles);          % Laser scan estimate from x_initial_estimate
    x_new_sm=zeros(3, num_particles);                     % New position after scan matching (ICP)
    x_sampled=zeros(3,num_samples, num_particles);        % Sample position around every particle position
    p_xt=zeros(num_samples, num_particles);               % Density function for odometry (MotionModelVelocity)
    p_zt=zeros(num_samples, num_particles);               % Density function for scans (Scan Matching)
    tau_xt=zeros(num_samples, num_particles);             % Optimal Proposal Distribution
    normalizer=zeros(1,num_particles);                    % Sample normalizer
    mu=zeros(3, num_particles);                           % Mean for Gaussian distribution from all particle samples 
    cov=zeros(3, 3, num_particles);                       % Covariance for Gaussian distribution from all particle samples 
    weight=zeros(size(weight_past, 1),1);                 % Particle weights
    x_final=zeros(3, num_particles);                      % Final estimate particle position after RBPF
    
%% RaoBlack Particle Filter for every particle
    for i = 1:num_particles
        %Estimate particle position from xt-1 after Ut
        con_cmd_noise = con_cmd+sqrt(PF.CovOdom)*randn(length(con_cmd),1);
        x_initial_estimate(:,i) = MotionCommandModel(x_past(:,i),con_cmd_noise,ts);
        x_initial_estimate(3,i) = NormalizeAngle(x_initial_estimate(3,i));
        
        %Estimation of laser given an particle position
        laser_estimate(:,:,i) = MeasurementEstimate(x_initial_estimate(:,i), PF.Probability(:,:,i),...
            LRF.AngleResolution, LRF.FOV, LRF.MaxDistance, PF.UsableArea);
             
        % Scan pre-processing prior to ICP algortithm, Consider Usable Area only
        true_scan = laser_data;
        estimate_scan = laser_estimate(:,:,i);
        
% DEBUG---------------------------------------------------------        
%         true_scan_debug = LocalPolar2World(x_true, true_scan);
%         estimate_scan_debug = LocalPolar2World(x_true, estimate_scan);
%             figure(fig_handle)
%             subplot(1,3,3)
%             caxis([0, 1]);
%             colormap gray;
%             colormap(flipud(colormap));
%             hold on;
%             imagesc(flipud(PF.Probability(:,:,i)));
%             % Plot scans
%             plot(true_scan_debug(1,:), true_scan_debug(2,:), '.g')
%             plot(estimate_scan_debug(1,:), estimate_scan_debug(2,:), '.b')
%----------------------------------------------------------------------------       
        
        true_scan(true_scan > PF.UsableArea) = nan;
        estimate_scan(estimate_scan > PF.UsableArea) = nan;

        % Chose usable data in both scans
        find_true_scan = isnan(true_scan(2,:));
        find_estimate_scan = isnan(estimate_scan(2,:));
        del_cols = or(find_true_scan, find_estimate_scan);
        true_scan(:,(del_cols(1,:)==1))=[];
        estimate_scan(:,(del_cols(1,:)==1))=[];

        % ICP returns the translation and rotation matrix and result data
        true_scan_pos=LocalPolar2World(x_true, true_scan);
        estimate_scan_pos=LocalPolar2World(x_initial_estimate(1:3,i), estimate_scan);
        failure=false;  
        try
            [TR, TT, modified_scan_pos] = ICP(true_scan_pos, estimate_scan_pos);       
        catch
            failure=true;
            disp('ICP Failed')
        end 
        if abs(TT(1))>50  || abs(TT(2))>50
            failure=true; % Check if ICP failed
        end
        % If ICP failed used last known position estimate
        if failure
            disp('ICP failed')
            x_final(:,i)=x_initial_estimate(:,i);  
            weight(i)=weight_past(i); % Keep particle weights unchanged                       
        else
            % Modify the pose according to ICP
            theta_new_sm = x_initial_estimate(3,i) + atan2(TR(2,1),TR(1,1));
            theta_new_sm = NormalizeAngle(theta_new_sm);
            position_new_sym = TR*x_initial_estimate(1:2,i) + TT;
            x_new_sm(:,i) = [position_new_sym; theta_new_sm];      
%% Plot grid-map, all positions and all scans
            figure(fig_handle)
            subplot(1,3,3)
            caxis([0, 1]);
            colormap gray;
            colormap(flipud(colormap));
            hold on;
            imagesc(flipud(PF.Probability(:,:,i)));
            % Plot scans
            plot(true_scan_pos(1,:), true_scan_pos(2,:), '.g')
            plot(estimate_scan_pos(1,:), estimate_scan_pos(2,:), '.b')
            plot(modified_scan_pos(1,:), modified_scan_pos(2,:), '.r')
            % Plot positions
            plot(x_true(1), x_true(2), 'dk')
            plot(x_initial_estimate(1,i), x_initial_estimate(2, i), 'db')
            plot(x_new_sm(1,i), x_new_sm(2,i), 'dr')
            axis([0 300 0 300])
            axis square
            hold off;
            pause(0.0001)
%% Randomize samples around good estimate (x_new_sm) of particle  
            x_sampled(:,:, i) = normrnd([x_new_sm(1,i).*ones(1, num_samples); ...
                                         x_new_sm(2,i).*ones(1, num_samples); ...
                                         x_new_sm(3,i).*ones(1, num_samples)], ...
                                         [0.1*ones(2, num_samples); 0.01*ones(1, num_samples)],...
                                         [3 num_samples]); % Choose Covarience Appropriately
            for s=1:num_samples           
               % Calculate velocity probability for that sample 
                p_xt(s,i)= MotionModelVelocity(x_sampled(:,s,i), con_cmd, x_past(:,i), a, ts);
               % Calculate scan matching probability for that sample 
                laser_true=laser_data;
                laser_est = MeasurementEstimate(x_sampled(:,s,i),PF.Probability(:,:,i),LRF.AngleResolution,LRF.FOV,LRF.MaxDistance,PF.UsableArea);
                laser_true(laser_true>PF.UsableArea)=nan;
                laser_est(laser_est>PF.UsableArea)=nan;
                find_laser_true=isnan(laser_true(2,:));
                find_laser_est=isnan(laser_est(2,:));
                del_cols = or(find_laser_true, find_laser_est);
                laser_est(:,(del_cols(1,:)==1))=[];
                laser_true(:,(del_cols(1,:)==1))=[];
                p_zt(s,i) = MeasurementModel(laser_true,laser_est,LRF.Sigma);
               % Calculate optimal proposal distribution
                tau_xt(s,i)=abs(p_xt(s,i)*p_zt(s,i));
            end
           % Calculate sample normalizer
            normalizer(i)=sum(tau_xt(:,i));
           % Calculate mu for gaussian distribution from samples
            mu(1, i)=(1/normalizer(i))*(x_sampled(1,:,i)*tau_xt(:,i));   
            mu(2, i)=(1/normalizer(i))*(x_sampled(2,:,i)*tau_xt(:,i));
            mu(3, i)=(1/normalizer(i))*(x_sampled(3,:,i)*tau_xt(:,i));
            mu(3, i)=NormalizeAngle(mu(3, i));
           % Calculate covariance for gaussian distribution from samples
            cov_sum=zeros(3,3);
            for s=1:num_samples
                cov_sum=cov_sum+(x_sampled(1:3,s,i)-mu(1:3,i))*(x_sampled(1:3,s,i)-mu(1:3,i))'*tau_xt(s,i);
            end
            cov(:,:,i)=(1/normalizer(i))*cov_sum;
           % Calculate new weights
            weight(i)=weight_past(i)*normalizer(i);
           % Use the Mean as final position of particle
            x_final(:,i)=mu(1:3,i);
        end
       % Make sure no errors have occured.
        if isnan(x_final(1,i)) || isnan(x_final(2,i))
                x_final(:,i)=x_initial_estimate(:,i);
                fprintf('Step %d, Particle %d X_state = NaN. Keep Unchanged. \n', iter, i);
        end
       % Update Grid Map for that particle
        [PF.Probability(:,:,i), PF.Log(:,:,i)] = GridMapping(PF.Log(:,:,i), x_final(:,i), laser_data, PF.UsableArea, PF.Log0, PF.GridSize, LRF);
    end
 
%% Resampling step
    weight=weight/norm(weight,1);            % Normalize weights
    Neff = 1/sum(weight.^2);                 % Calculate Neff for possible resampling
    resample_percentage = 0.5;               % 0.5 is the best resample percentage proposed by Cyrill Stachniss
    Nt = resample_percentage*num_particles;
    if Neff < Nt
%         disp('Resampling');
        [x_final(:,:), weight, idx] = Resample(x_final(:,:), weight);
        P_save_tmp=PF.Probability(:,:,:);
        L_save_tmp=PF.Log(:,:,:);
        for s=1:num_particles
           PF.Probability(:,:,s)=P_save_tmp(:,:,idx(s));
           PF.Log(:,:,s)=L_save_tmp(:,:,idx(s));
        end
    end
    PF.Weight(:,iter) = weight;                      % Save particle weights (might have changed during resampling)  
    PF.Particles(:,:,iter) = x_final(:,:);           % Save particle updates
    [~, max_index] = max(PF.Weight(:,iter));
    PF.Best(iter) = max_index;                       % Logging the Best Particle
%% Calculate estimate position given all weighted particles
    x_state = zeros(3,1);
    for i=1:num_particles
            x_state = x_state +(weight(i).*x_final(:,i));
    end
    x_state(3) = NormalizeAngle(x_state(3));
    PF.Xstate(:, iter)= x_state;
end