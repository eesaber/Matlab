function [p] = Stat(A)
    % STAT function is used to scale the image.
    % USAGE:
    % [p] = STAT(A), A is the image matrix and p is the 1-d array, which 
    % can be recovered to image size by using RESHAPE(p, SIZE(A)). The
    % return value is unsigned 8-bit integer.
    % ALGORITHM:
    % 1. Turn matrix A to dB unit.
    % 2. Determines the mean and standard deviation of whole image.
    % 3. Map range [μ-2σ, μ+2σ] to [0 255], set the value which is smaller
    %    than μ-2σ to μ-2σ; and set the value which is larger than μ+2σ to
    %    μ+2σ.
    %....................................................................%
    
    p = reshape(A,[1,numel(A)]);
    p_sigma = std(p(~isinf(p)),'omitnan');
    p_mean = mean(p(~isinf(p)),'omitnan');
    p(p > p_mean+2*p_sigma) = p_mean+2*p_sigma;
    p(p < p_mean-2*p_sigma) = p_mean-2*p_sigma;
    p = uint8((p - (p_mean-2*p_sigma))*255/(4*p_sigma));
end