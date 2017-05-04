function [ omega ] = omega_fun(scene)
    if sum(sum(scene)) == 0
        omega = 0;
    elseif scene(2,2) == 1
        omega = (sum(sum(scene)) - 1) / 8;
    else
        omega = (sum(sum(scene))) / 8;
    end
end