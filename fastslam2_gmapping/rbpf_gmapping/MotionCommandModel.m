function [x_new] = MotionCommandModel(x,u,ts)
%% Calculates the end position given an initial position and a movemnt command
    if u(2)==0
        u(2)=1e-10;     % If rotational value of command is exactly zero it gives errors
    end
    x_new = x + [-(u(1)/u(2))*sin(x(3))+(u(1)/u(2))*sin(x(3)+u(2)*ts); ...
                  (u(1)/u(2))*cos(x(3))-(u(1)/u(2))*cos(x(3)+u(2)*ts); ...
                   u(2)*ts];
    x_new(3) = NormalizeAngle(x_new(3));
end