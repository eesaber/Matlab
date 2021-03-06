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
    %   The implementation is copied from the python lib "Keras", 
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
    %% Initialization of parameters
    %rng('shuffle')
    rng(87)
    im_h = size(im,1); % original width 
    im_w = size(im,2); % original height
    img_h = 256; % shrinked height size
    img_w = 256; % shrinked width size
    shear_range = p_.shear_range;
    rotation_range = p_.rotation_range;
    zoom_range = p_.zoom_range;
    horizontal_flip = p_.horizontal_flip;
    vertical_flip = p_.vertical_flip;
    OUTPUT_PATH = [obj.OUTPUT_PATH, '/aug/'];
    enhence = 40;
    x_train = single(zeros(enhence+p_.data_size, img_h, img_w, size(im,3)));
    y_train = single(zeros(enhence+p_.data_size, img_h, img_w));
    %% Data augmentation
    % [1 shy 0; shx 1 0; 0 0 1]
    % where shx and shy specifies the shear factor along the x axis and y
    % axis, respectively.
    for it = 1 : p_.data_size
        theta = rand(1)*rotation_range;
        shy = rand(1)*shear_range;
        shx = rand(1)*shear_range;
        tform = affine2d(...
            [1, shy, 0;...
            shx, 1, 0;...
            0, 0, 1]);
        %cut = randi([1 im_w]);
        %temp_im = cat(2, im(:,cut+1:end,:), im(:,1:cut,:));
        %temp_label = [label(:,cut+1:end), label(:,1:cut)]; 
        temp_im = im;
        temp_label = label;
        zoom_w = max(img_w, int32((zoom_range*rand(1))*im_w));
        zoom_h = max(img_h, int32((zoom_range*rand(1))*im_h));
        col_ = randi([1, im_w-zoom_w+1],1);
        row_ = randi([1, im_h-zoom_h+1],1);
        temp_im = temp_im(row_: row_+zoom_h-1, col_: col_+zoom_w-1, :);
        temp_label = temp_label(row_: row_+zoom_h-1, col_: col_+zoom_w-1);
        
        temp_im = permute(imresize(imwarp(temp_im,tform), [img_h, img_w]),[1,2,3]);
        temp_im = imrotate(temp_im, theta, 'bilinear','crop');
        temp_label = permute(imresize(imwarp(temp_label,tform), [img_h, img_w]),[1,2]);
        temp_label = imrotate(temp_label, theta, 'nearest','crop');
        if vertical_flip && logical(randi([0 1]))
            temp_im = flipud(temp_im);
            temp_label = flipud(temp_label);
        end
        if horizontal_flip && logical(randi([0 1]))
            temp_im = fliplr(temp_im);
            temp_label = fliplr(temp_label);
        end
        x_train(it,:,:,:) = temp_im;
        y_train(it,:,:) = temp_label;
        
        imwrite(temp_label,[OUTPUT_PATH, num2str(it), '.png'])
        %{
        figure
        subplot(2,1,1)
        image(cat(3,temp_im, ones(img_h, img_w)))
        subplot(2,1,2)
        imagesc(temp_label)
        drawnow
        %}
    end
    
    % Enhence on the training on the deformed FYI
    %im(1:264, 1:826)
    for it = p_.data_size+1 : p_.data_size+enhence        
        theta = rand(1)*rotation_range;
        shy = rand(1)*4*shear_range;
        shx = rand(1)*4*shear_range;
        tform = affine2d(...
            [1, shy, 0;...
            shx, 1, 0;...
            0, 0, 1]);
        
        temp_im = im(1:264, 1:826,:);
        temp_label = label(1:264, 1:826);
        temp_im = permute(imresize(imwarp(temp_im,tform), [img_h, img_w]),[1,2,3]); % shear 
        temp_im = imrotate(temp_im, theta, 'bilinear','crop'); % rotate
        temp_label = permute(imresize(imwarp(temp_label,tform), [img_h, img_w]),[1,2]);
        temp_label = imrotate(temp_label, theta, 'nearest','crop');
        if vertical_flip && logical(randi([0 1])) % vertical flip
            temp_im = flipud(temp_im);
            temp_label = flipud(temp_label);
        end
        if horizontal_flip && logical(randi([0 1])) % horizontal flip
            temp_im = fliplr(temp_im);
            temp_label = fliplr(temp_label);
        end
        x_train(it,:,:,:) = temp_im;
        y_train(it,:,:) = temp_label;
        %imwrite(temp_label,[OUTPUT_PATH, num2str(it), '.png'])
    end
    
    %% Save results
    save([obj.OUTPUT_PATH, '/x_train_', p_.filename], 'x_train')
    save([obj.OUTPUT_PATH, '/y_train_', p_.filename], 'y_train')
    % Output results
    if nargout>=1, varargout{1} = x_train; end
    if nargout>=2, varargout{2} = y_train; end
    
end