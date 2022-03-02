function [value,isterminal,direction] = Devents(t,y)
    value = (y(1)^2 - 1) * (y(2)^2 - 1);
    isterminal = 1;
    direction = 0;
end