function p_xt = MotionModelVelocity(x_sample, con_cmd, x_past, a, ts) 
%% Calculates the probability desity function of a specific position given an initial position and a movement command
    x_old=x_past(1);
    y_old=x_past(2);
    theta_old=NormalizeAngle(x_past(3));

    x_new=x_sample(1);
    y_new=x_sample(2);
    theta_new=NormalizeAngle(x_sample(3));

    v=con_cmd(1);
    w=rad2deg(con_cmd(2));
    
    % Motion Model Velocity Algorithm in Probabilistic Robotics Book
    mu= 0.5*(((x_old-x_new)*cos(theta_old)+(y_old-y_new)*sin(theta_old))/...
        ((y_old-y_new)*cos(theta_old)-(x_old-x_new)*sin(theta_old)));
    
    if mu==Inf
        mu=-abs(x_old-x_new);
    end
    if mu==-Inf
        mu=abs(x_old-x_new);
    end

    x_star=((x_old+x_new)/2)+mu*(y_old-y_new);
    y_star=((y_old+y_new)/2)+mu*(x_new-x_old);
    r_star=sqrt((x_old-x_star)^2+(y_old-y_star)^2);
    theta_diff=atan2((y_new-y_star),(x_new-x_star))-atan2((y_old-y_star),(x_old-x_star));
    theta_diff=NormalizeAngle(theta_diff);

    v_new=(theta_diff/ts)*r_star;
    w_new=theta_diff/ts;
    % theta_2=theta_1+w_new;
    ProbNormDistribution = @(a,b)(1/sqrt(2*pi*b))*exp(-0.5*(a^2/b));
    p_xt1=ProbNormDistribution(v-v_new,a(1)*v^2+a(2)*w^2);
    p_xt2=ProbNormDistribution(w-w_new,a(3)*v^2+a(4)*w^2);
    gamma=((theta_new-theta_old)/ts)-w_new;
    p_xt3=ProbNormDistribution(gamma,a(5)*v^2+a(6)*w^2);      
%   Do Not Consider The Last Rotation
    p_xt=p_xt1*p_xt2;
end