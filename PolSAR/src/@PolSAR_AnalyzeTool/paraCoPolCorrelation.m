function co_pol = paraCoPolCorrelation(obj, x, y)
    % PARACOPOLCORRELATION Co-polarization correlation
    %
    % Syntax:  PARACOPOLCORRELATION(x,y)
    %
    % Inputs:
	%    x  - row number of boxcar filter with default value 3.
	%    y  - column number of boxcar filter with default value 3.
    %
    % Other m-files required: plot_para
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    switch nargin
        case 3
            mask = ones(x,y)/(x*y);
        case 2
            mask = ones(x,3)/(x*3);
        otherwise
            mask = ones(3,3)/(3*3);
    end
    if numel(obj.vv_vv)==0 || numel(obj.hh_hh)==0
    error('vv_vv can not be empty')
    end

    disp('S_{hh} S_{vv}^* / \sqrt{σ_vv σ_hh}')
    f_name = 'para_copolcorrelation';
    co_pol = (angle(conv2(obj.hh_vv./sqrt(obj.hh_hh.*obj.vv_vv), mask, 'same')));
    %co_pol = abs(conv2(obj.hh_vv, mask, 'same')./ ...
    %    sqrt(conv2(obj.hh_hh, mask, 'same').*conv2(obj.vv_vv, mask, 'same')));
    %{
    %co_pol = 180/pi*angle(conv2(obj.hh_vv./sqrt(obj.hh_hh.*obj.vv_vv), mask, 'same'));
    %co_pol = 10*log10(obj.hh_hh);
    co_pol = rescale(co_pol, 0, 1);
    co_pol = histeq(co_pol);
    co_pol = adapthisteq(co_pol);
    
    co_pol = entropyfilt(co_pol);
    co_pol = imgaussfilt(co_pol, 8);
    
    co_pol = rescale(co_pol, 0, 1);
    co_pol = histeq(co_pol);
    %co_pol = imsharpen(co_pol);
    %}
    figure
    imagesc(rescale((co_pol)))
    obj.plotSetting([0 1])
    plot_para('Maximize',true, 'Ratio', [4 3 1], 'Filename',[obj.OUTPUT_PATH '/' f_name])
    
end