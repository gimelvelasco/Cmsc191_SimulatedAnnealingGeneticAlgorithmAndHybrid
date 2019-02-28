function F_obj_ret = OBJFUNC(x,tpl,test_func)   %Objective: The Fitness Function must be very close if not equal to 0
    if test_func == 1       %DeJong
        F_obj_ret = 0;
        for i=1:tpl
           F_obj_ret = F_obj_ret + x(i)^2;
        end
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 2   %Axis Parallel Hyper-ellipsoid
        F_obj_ret = 0;
        for i=1:tpl
           F_obj_ret = F_obj_ret + i*x(i)^2;
        end
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 3   %Rotated Hyper-ellipsoid
        F_obj_ret = 0;
        for i=1:tpl
            for j=1:i
                F_obj_ret = F_obj_ret + x(j)^2;
            end
        end
        F_obj_ret = abs(F_obj_ret);
    elseif test_func == 4   %Rastrigin
        F_obj_ret = 0;
        for i=1:tpl
           F_obj_ret = F_obj_ret + x(i)^2 - 10*cos(2*pi*x(i));
        end
        F_obj_ret = F_obj_ret + 10*tpl;
        F_obj_ret = abs(F_obj_ret);
    else   %Ackley
        a = 20;
        b = 0.2;
        c = 2*pi;
        F_obj_ret1 = 0;
        F_obj_ret2 = 0;
        for i=1:tpl
            F_obj_ret1 = F_obj_ret1 + x(i)^2;
            F_obj_ret2 = F_obj_ret2 + (cos(c*x(i)));
        end
        F_obj_ret1 = (-a)*exp((-b)*sqrt((1/tpl)*F_obj_ret1));
        F_obj_ret2 = exp((1/tpl)*F_obj_ret2);
        F_obj_ret1 = F_obj_ret1 - F_obj_ret2;
        F_obj_ret2 = a + exp(1);
        F_obj_ret = F_obj_ret1 + F_obj_ret2;
        F_obj_ret = abs(F_obj_ret);
    end
end