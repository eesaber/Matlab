function varargout = myDataAugmentation(obj, im, label, varargin)
    % MYDATAUGMENTATION implements data augmentation for sea-ice
    % classification.
    %
    % Syntax: 
    % * MYDATAUGMENTATION(im)
    %
    % Inputs:
    %    * im - (n,m,d) array
    %    
    % Example: 
    %    1. GETPOLANGLE([0 90])
    %       Set the angle range = [0 90].
    % Reference 
    %   The implementation is copFrom the python lib "Keras", 
    %   keras.preprocessing.image.ImageDataGenerator(
    %   featurewise_center = False,
    %   samplewise_center=False, 
    %   featurewise_std_normalization=False,
    %   samplewise_std_normalization=False,
    %   zca_whitening=False,
    %   zca_epsilon=1e-06, 
    %   rotation_range=0,
    %   width_shift_range=0.0,
    %   height_shift_range=0.0,
    %   brightness_range=None,
    %   shear_range=0.0,
    %   zoom_range=0.0,
    %   channel_shift_range=0.0,
    %   fill_mode='nearest',
    %   cval=0.0,
    %   horizontal_flip=False,
    %   vertical_flip=False, 
    %   rescale=None)
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    
    %% Parse input and output arguments
    p_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x, {'numeric'}, {});
    validationFcn_2_ = @(x) validateattributes(x, {'logical'}, {});
    validationFcn_3_ = @(x) validateattributes(x, {'char'}, {'nonempty'});
    addParameter(p_, 'data_size', 32, validationFcn_1_);
	addParameter(p_, 'shear_range', 0, validationFcn_1_);
    addParameter(p_, 'rotation_range', 0, validationFcn_1_);
    addParameter(p_, 'zoom_range', 0, validationFcn_1_);
    addParameter(p_, 'horizontal_flip', false, validationFcn_2_);
    addParameter(p_, 'vertical_flip', false, validationFcn_2_);
    addParameter(p_, 'filename', 'augmentation_data', validationFcn_3_);
	parse(p_,varargin{:})
    p_ = p_.Results;
    
    % parse output arguments
    minArgs=0;
    maxArgs=2;
    nargoutchk(minArgs,maxArgs)
    %%
    
    close all
    im = imread('peppers.png');
    
    im_h = size(im,1);
    im_w = size(im,2);
    img_h = 496;
    img_w = 496;
    label = zeros(img_h, img_w);
    %%
    shear_range = p_.shear_range;
    rotation_range = p_.rotation_range;
    zoom_range = p_.zoom_range;
    x_train = single(zeros(p_.data_size, img_h, img_w, size(im,3)));
    y_train = single(zeros(p_.data_size, img_h, img_w));
    % [1 shy 0; shx 1 0; 0 0 1]
    % where shx and shy specifies the shear factor along the x axis and y
    % axis, respectively.
   
    if p_.horizontal_flip
        horizontal_flip = logical(randi([0,1], p_.data_size,1));
    else
        horizontal_flip = false(p_.data_size,1);
    end
    if p_.vertical_flip
        vertical_flip = logical(randi([0,1], p_.data_size,1));
    else
        vertical_flip = false(p_.data_size,1);
    end
    
    for it = 1 : p_.data_size
        theta = rand(1)*rotation_range;
        shy = rand(1)*shear_range;
        shx = rand(1)*shear_range;
        tform = affine2d(...
            [1, shy, 0;...
            shx, 1, 0;...
            0, 0, 1]);
        
        zoom = 1-zoom_range*rand(1);
        zoom_w = int32(zoom*im_w);
        zoom_h = int32(zoom*im_h);
        col_ = randi([1, im_w-zoom_w+1],1);
        row_ = randi([1, im_h-zoom_h+1],1);
        temp_im = im(row_: row_+zoom_h-1, col_: col_+zoom_w-1, :);
        temp_label = label(row_: row_+zoom_h-1, col_: col_+zoom_w-1, :);
        
        temp_im = permute(imresize(imwarp(temp_im,tform), [img_h, img_w]),[1,2,3]);
        temp_im = imrotate(temp_im, theta, 'bilinear','crop');
        temp_label = permute(imresize(imwarp(temp_label,tform), [img_h, img_w]),[1,2]);
        temp_label = imrotate(temp_label, theta, 'nearest','crop');
        if vertical_flip(it) 
            temp_im = flipud(temp_im);
            temp_label = flipud(temp_label);
        end
        if horizontal_flip(it) 
            temp_im = fliplr(temp_im);
            temp_label = fliplr(temp_label);
        end
        x_train(it,:,:,:) = temp_im;
        y_train(it,:,:) = temp_label;
        figure
        image(temp_im)
        drawnow
        
    end
    save([obj.OUTPUT_PATH, '/x_train_', p_.filename], 'x_train')
    save([obj.OUTPUT_PATH, '/y_train_', p_.filename], 'y_train')
    % Output results
    if nargout>=1, varargout{1} = x_train; end
    if nargout>=2, varargout{2} = y_train; end
    
end