function [x_world] = LocalPolar2World(x_robot, polar)
 %% Cnverts polar coordinats [alpha;r] to world [x_world]
    alpha = polar(1,:);
    r = polar(2,:);
    x = zeros(1,size(polar,2));
    y = zeros(1,size(polar,2));

    for i = 1:size(polar,2)
        x(i) = r(i)*cos(alpha(i));
        y(i) = r(i)*sin(alpha(i));
    end
    local_cart = [x;y];
    x_world = zeros(size(local_cart));
    for i = 1:size(local_cart,2)
        x_world(1,i) = cos(x_robot(3))*(local_cart(1,i)) - sin(x_robot(3))*(local_cart(2,i)) + x_robot(1);
        x_world(2,i) = cos(x_robot(3))*(local_cart(2,i)) + sin(x_robot(3))*(local_cart(1,i)) + x_robot(2);
    end
end