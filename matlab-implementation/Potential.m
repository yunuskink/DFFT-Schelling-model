function [V_p]= Potential(pi,P)
% Define Spatial Utility
    if P(pi,3)==1 % red
        x=P(pi,1);
        y=P(pi,2);
        
        V_p=-1.3*y; % Adjust This Part
    else % blue
        x=P(pi,1);
        y=P(pi,2);
        
        V_p=-1.3*x; % Adjust This Part
    end
end