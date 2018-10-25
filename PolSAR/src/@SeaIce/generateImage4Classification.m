function im = generateImage4Classification(obj, texture)

    switch nargin
        case 1
            % Produce texture paramters
            temp = 10*log10(obj.vv_vv);
            temp(10*log10(obj.vv_vv)<-25) = -25;
            temp(10*log10(obj.vv_vv)>-5) = -5;
            texture = myGLCM(rescale(temp),21,41);
            clear temp
    end
    
    %
    obj.logCumulant();
    if isempty(obj.T_11), obj.cov2coh(); end

    pol_ang = obj.getPolAngle([-30 30]);
    poolang = obj.myMaxPooling(abs(pol_ang),9,15);
    poolang = imgaussfilt(imresize(poolang,obj.IMAGE_SIZE,'nearest'),3);

    %% input image tuning  
    input_opt = 4;
    switch input_opt
        case 1
            im = reshape(cat(3, rescale(obj.kai_1), rescale(obj.kai_2),...
            rescale(obj.kai_3), rescale(obj.kai_4), 0.3*rescale(abs(poolang))),[],5);
        case 2
            im = reshape(cat(3, rescale(obj.kai_1), rescale(obj.kai_2),...
            rescale(obj.kai_3), 0.3*rescale(poolang)),[],4);
        case 3
            im = reshape(cat(3, rescale(obj.kai_1), rescale(obj.kai_2),...
            rescale(obj.kai_3), rescale(obj.kai_4)),[],4);
        case 4 
            im = reshape(cat(3, rescale(obj.kai_1), rescale(obj.kai_2),...
                rescale(texture(:,:,1))),[],3);
        case 5
            im = reshape(rgb2lab(cat(3, rescale(obj.kai_1), rescale(obj.kai_2),...
                rescale(texture(:,:,1)))),[],3);
            im = im(:,2:3);
        otherwise
            im = reshape(cat(3, obj.kai_1, obj.kai_2,obj.kai_3, obj.kai_4),[],4);
    end
end