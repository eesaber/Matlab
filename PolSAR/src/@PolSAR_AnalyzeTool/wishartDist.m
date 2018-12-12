function dist_ = wishartDist(obj, x, y)
    % WISHARTDIST implement the revised wishart distance.
    %
    % Inputs:
    %   * x is a (3,3,n) ndarray 
    %   * y is a (3,3,m) ndarray 
    %
    % Definition:
    %   * D(a, b) = log(|b|) + Tr(b^{-1}*a)
    %     where base of log is 10, a is observation and b is cluster centroid.
    %   * D(a, b) = log(|b|) - log(|a|) + Tr(b^{-1}*a)
    % Reference:
    % [1] Unsupervised classification using polarimetric decomposition
    %     and the complex Wishart classifier.
    % [2] Unsupervised classification of polarimetric synthetic aperture
    %     radar images using fuzzy clustering and EM clustering
    %------------- <<<<< >>>>>--------------
    switch nargin
        case 2
            dist_ = pairdist(x);
        case 3
            dist_ = pairdist2(x,y);
        otherwise 
            error('At least one input arguments required.')
    end
end

function dist_ = pairdist(x)
    dist_ = single(zeros((size(x,3)-1)/2*size(x,3),1));
    counter = 1;
    size_x = size(x,3);
    for it_i = 1 : size_x
        ln_det_i = log10(real(det(x(:,:,it_i))));
        inv_i = inv(x(:,:,it_i));
        for it_j = it_1+1 : size_x
            ln_det_j = log10(real(det(x(:,:,it_j))));
            dist_(counter) = ln_det_i-ln_det_j+...
                real(trace(inv_i*x(:,:,it_j)));
            counter = counter+1;
        end
    end
end

function dist_ = pairdist2(x, y)
    % Maintain total number of y is less than total number of x
    flag_exg = false;
    if size(y,3) > size(x,3)
        flag_exg = true;
        temp = y;
        y = x;
        x = temp;
        clear temp
    end
    % dist_(i,j) is the distance between x(:,:,i) and y(:,:,j)
    dist_ = single(zeros(size(x,3), size(y,3)));
    % calculate distance
    parfor it_y = 1 : size(y,3)
        ln_det_y = log10(real(det(y(:,:,it_y))));
        inv_y = inv(y(:,:,it_y));
        tmp = single(zeros(size(x,3),1));
        for it_x = 1 : size(x,3)
            ln_det_x = log10(real(det(x(:,:,it_x))));
            tmp(it_x) = ln_det_y-ln_det_x+...
                real(trace(inv_y*x(:,:,it_x)));
        end
        dist_(:,it_y) = tmp;
    end
    if flag_exg 
        dist_ = dist_.';
    end
end