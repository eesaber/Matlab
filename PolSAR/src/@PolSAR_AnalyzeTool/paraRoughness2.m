function paraRoughness2(obj, x, y)
    % PARAROUGHNESS2 Circular basis co-polarization coefficient (T_22+T_33)/(T_22 -T_33)
    % The ratio relates to surface roughness and is independent to 
    % the permittivity of surface for ks < 1 where k is radio wavenumber and s is 
    % the rms height of surface [1] [2] [3].
    %
    % Syntax: PARAROUGHNESS2(x,y)
    %
    % Inputs:
	%    x  - row number of boxcar filter with default value 3.
	%    y  - column number of boxcar filter with default value 3.
    %
    % Other m-files required: plot_para
    % Subfunctions: none
    % MAT-files required: none
    %
    % Reference 
    % [1] The effect of surface roughness on multifrequency polarimetric SAR data
    % [2] Surface roughness and slope mesaruements using polarimetric SAR data
    % [3] Inversion of surface parameters from polarimetric SAR
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    disp('Surface roughness by \frac{T_22 + T_33}{T_22 - T_33}}')
    switch nargin
    case 3
        mask = ones(x,y)/(x*y);
    case 2
        mask = ones(x,3)/(x*3);
    otherwise
        mask = ones(3,3)/(3*3);
    end
    if numel(obj.T_11)==0 || numel(obj.T_22)==0 || numel(obj.T_12)==0
        error('T_11 can not be empty')
    end
    f_name = 'para_roughness1';
    figure
    imagesc(conv2(abs(obj.T_12)./sqrt(obj.T_11.* obj.T_22), mask, 'same'))
    obj.plotSetting([0 1])
    plot_para('Maximize',true,'Ratio', [4 3 1], 'Filename',[obj.OUTPUT_PATH '/' f_name])
end