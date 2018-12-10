function [im,texture] = generateImage4Classification(obj, input_vector, varargin)
    % GENERATIONIMAGE4CLASSIFICATION generate the input image for sea-ice classification.
    %
    % syntax:
    % GENERATIONIMAGE4CLASSIFICATION(intput_vecor, name, value)
    % Specifies function properties using one or more Name,Value pair arguments. 
    %
    % Inputs:
    % * input_vector: should be an interger form 1 to 4.
    %
    % * Name-Value Pair Arguments:
    %    ('texture', texture) - The pixel above SNR thershold is marked. If s is ture, 
    %                function will show the figure.
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    
    p_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'numeric'},{'nonempty'});
    addParameter(p_,'texture',[],validationFcn_1_);
    parse(p_,varargin{:})
    p_ = p_.Results;
    
    if isempty(p_.texture) && (input_vector==1 || input_vector==2)
        temp = 10*log10(obj.vv_vv);
        temp(10*log10(obj.vv_vv)<-25) = -25;
        temp(10*log10(obj.vv_vv)>-5) = -5;
        texture = myGLCM(rescale(temp),21,41);
        %texture = myGLCM(rescale(temp),11,21);
        clear temp
    else
        texture = p_.texture;
    end
    
    if input_vector~=4 && (isempty(obj.kai_1) || isempty(obj.kai_2) || isempty(obj.kai_3) || isempty(obj.kai_4))
        obj.logCumulant();
    end
    if isempty(obj.H) && input_vector==2
        [~,~] = obj.eigenDecomposition(1,'eigen_9x9avg','SaveResults', false);
    end
    %% input vector generation
    im = reshape(cat(3, rescale(obj.kai_1,'InputMax',-5,'InputMin',-13),...
                rescale(obj.kai_2,'InputMax',7,'InputMin',0)),[],2);
    switch input_vector
        case 1
            im = [im, reshape(rescale(texture(:,:,1)),[],1)];
        case 2
            co_pol = 10*log10(obj.vv_vv./obj.hh_hh);
            co_pol(co_pol>3) = 3;
            co_pol(co_pol<-3) = -3;
            co_pol = rescale(co_pol);    
            im = [im, reshape(rescale(texture(:,:,1)),[],1), obj.H(:),...
                co_pol(:)];
            clear co_pol
        case 3
            1;
        case 4
            bai = 100;
            im = bai*reshape(cat(3, obj.hh_hh, obj.hv_hv, obj.vv_vv,...
                real(obj.hh_hv), imag(obj.hh_hv),...
                real(obj.hh_vv), imag(obj.hh_vv),...
                real(obj.hv_vv), imag(obj.hv_vv)),[],9);
        case 5
            %{
            ang = 180/pi*obj.getPolAngle([-30 30]);
            ang(ang>30) = 30;
            ang(ang<-30) = -30;
            ang = rescale(ang);
            im = [im, ang(:)];
            %}
            temp = 10*log10(obj.hv_hv);
            temp(temp<-30) = -30;
            temp(temp>-15) = -15;
            temp = rescale(temp);
            im = [im, temp(:)];
    end
end