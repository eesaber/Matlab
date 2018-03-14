function [theta] = Find_angle(T)
% FIND_ANGLE implements the algorithm [1] that obtain the oreintatation angle 
% from covaraince matrix C
% 
% [1] "On the Estimation of Radar Polarization Orientation Shifts Induced by Terrain Slopes"
    [~,~,size_T] = size(T);
    theta = 1/4*(atand(-4*real(T(2,3,:))./2./(-T(2,2,:)+4*T(3,3,:)))+pi);
    theta(theta > pi/4) = theta(theta > pi/4) - pi/2;
end