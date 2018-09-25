%% Function Definition
function [dr] = beadWire(~,y)
    global OMEGA m g r;
    dr = zeros(2,1);
    
    dr(1) = y(2);
    dr(2) = (-g/r)*cos(y(1)) - OMEGA^2*cos(y(1))*sin(y(1));
end