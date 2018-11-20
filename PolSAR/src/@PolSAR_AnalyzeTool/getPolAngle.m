function theta = getPolAngle(obj, ang_range)
    % GETPOLANGLE implements the algorithm [1] that obtain the oreintatation angle 
    % from covaraince matrix C
    %
    % Syntax: 
    % * GETPOLANGLE(ang_range)
    %
    % Inputs:
    %    * ang_range - 1x2 array, Specify angle of limits as a 
    %      two-element vector of the form [xmin xmax], where xmax is 
    %      greater than xmin.
    %    
    % Example: 
    %    1. GETPOLANGLE([0 90])
    %       Set the angle range = [0 90].
    %
    % Reference :
    % [1] "On the Estimation of Radar Polarization Orientation Shifts 
    %      Induced by Terrain Slopes"   
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    
    % atan2(Y,X), returns values in the closed interval [-pi,pi]
    % atan(X), returns values in the closed interval [-pi/2,pi/2]
    mask = ones(2,4)/8;
    if isempty(obj.T_11), obj.cov2coh(); end    
    %theta = 1/4*(atan2(-4*real(obj.T_23), -obj.T_22 + obj.T_33)+pi);
    theta = 1/4*(atan2(-4*real(conv2(obj.T_23,mask,'same')), -conv2(obj.T_22,mask,'same') + conv2(obj.T_33,mask,'same'))+pi);
    theta(theta > pi/4) = theta(theta > pi/4) - pi/2;
    figure
        imagesc(theta/pi*180)
        obj.plotSetting(ang_range,'Colorbar_unit',"(deg)")
        plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH '/polAngle'])
end