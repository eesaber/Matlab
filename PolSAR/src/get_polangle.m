function kaku = get_polangle(T, plot_set, ang_range)
% FIND_ANGLE implements the algorithm [1] that obtain the oreintatation angle 
% from covaraince matrix C
% 
% [1] "On the Estimation of Radar Polarization Orientation Shifts Induced by Terrain Slopes"
    global im_size
    % atan2(Y,X), returns values in the closed interval [-pi,pi]
    % atan(X), returns values in the closed interval [-pi/2,pi/2]
    theta = squeeze(1/4*(atan2(-4*real(T(2,3,:)), -T(2,2,:)+T(3,3,:))+pi)).';
    theta(theta > pi/4) = theta(theta > pi/4) - pi/2;
    sig_angle = reshape(theta, im_size);
    kaku = sig_angle/pi*180;
    figure
        imagesc(sig_angle/pi*180)
        plot_set(ang_range,1)
        plot_para('Filename','angle_sig_1', 'Maximize',true)
end
