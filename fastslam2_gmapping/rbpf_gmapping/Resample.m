function [x_new, weight_new, idx] = Resample(x, weight)
%% Resampling function
    Ns = length(weight);  % Ns = number of particles
    replacement = true;
    idx = randsample(1:Ns, Ns, replacement, weight);
    x_new = x(:,idx);                    % extract new particles
    weight_new = repmat(1/Ns, 1, Ns);    % now all particles have the same weight
end