function [X] = HMRF(obj, im, labels, nColors)
    %% HMRF implement hidden Markov random field.

    if [size(im,1) size(im,2)] ~= size(labels)
        error('The dimension of labels and im are different')
    end
    EM_iter=10; % max num of iterations
    MAP_iter=10; % max num of iterations
    mag = single(rescale(histeq(rgb2gray(im), 2^16-1)));
    for n = 1 : nColors
        mu(n) = mean(mag(labels==n));
        sigma(n) = std(mag(labels==n));
    end
    Z = edge(mag, 'Canny', 0.75);
    figure
    imagesc(Z)
    colormap gray
    set(gca, 'Ydir','normal')
    disp('Continue?')
    pause(5)
    tic
    [X, ~, ~] = HMRF_EM(labels, mag , Z, mu, sigma,...
                   nColors, EM_iter, MAP_iter);
    toc
end