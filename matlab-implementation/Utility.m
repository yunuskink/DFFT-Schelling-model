function [u]=Utility(NO,NG,og)
% Define Social Utility
    if og==0
        u=0.4*NO;
    else
        u=0.4*NG;
    end
end