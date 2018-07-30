function [p] = Stat(A, bit)
    % STAT function is used to scale the image.
    % USAGE:
    % 1. [p] = STAT(A), A is the image matrix and p is the 1-d array, which 
    % can be recovered to image size by using RESHAPE(p, SIZE(A)). The
    % return value is unsigned 8-bit integer.
    % 2. [p] = STAT(A, bit), bit can be 'uint8', 'uint16' or 'uint32'. The
    % defalut data type is uint8.
    % ALGORITHM:
    % 1. Turn matrix A to dB unit.
    % 2. Determines the mean and standard deviation of whole image.
    % 3. Map range [μ-2σ, μ+2σ] to [0 255], set the value which is smaller
    %    than μ-2σ to μ-2σ; and set the value which is larger than μ+2σ to
    %    μ+2σ.
    %....................................................................%
    switch nargin
        case 1
            scale = 2^8-1;
            f = @(x) uint8(x);
        case 2
            if strcmp('uint16',bit)
                scale = 2^16-1;
                f = @(x) uint16(x);
            elseif strcmp('uint32',bit)
                scale = 2^32-1;
                f = @(x) uint32(x);
            else
                scale = 2^8-1;
                f = @(x) uint8(x);    
            end
    end
    
    p = reshape(A,[1,numel(A)]);
    p_sigma = std(p(~isinf(p)),'omitnan');
    p_mean = mean(p(~isinf(p)),'omitnan');
    p(p > p_mean+2*p_sigma) = p_mean+2*p_sigma;
    p(p < p_mean-2*p_sigma) = p_mean-2*p_sigma;
    p = f((p - (p_mean-2*p_sigma))*scale/(4*p_sigma));
end