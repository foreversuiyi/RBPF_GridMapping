function q = MeasurementModel(laser_true, laser_est, sigma)
%% Calculates the probability of one point in the real scan matching to that same point in the estimated scan
% How many beams are considered
model_num = 2;
model_step = length(laser_true)/model_num;
q = 1;
sigma = 100 * sigma;  % This Multiplication is Important in case of q = 0;
% Simply Using Gaussian Measurement Model in Probabilistic Robotics Book
for i=1:model_num
    step = floor(model_step*i);
    p_a = exp(-0.5*((laser_true(1,step) - laser_est(1,step))/(sigma(1)))^2)/(sqrt(2*pi)*sigma(1));
    p_r = exp(-0.5*((laser_true(2,step) - laser_est(2,step))/(sigma(2)))^2)/(sqrt(2*pi)*sigma(2));
    p = p_a*p_r;
    q= q*p;
end